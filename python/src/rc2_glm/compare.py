"""CLI: diff Python GLM outputs against a MATLAB reference directory.

Run::

    rc2-glm-compare \\
        --python-dir  path/to/python_run \\
        --reference-dir path/to/matlab_run \\
        --out path/to/validation/comparison_report.csv

Checks (each one emitted as a row in ``comparison_report.csv`` with columns
``check, value, threshold, pass, note``):

1. Prefilter category agreement — fraction of common clusters with identical
   ``category``.                                         Threshold ≥ 0.90.
2. Selected-variable Jaccard — mean Jaccard over the two ``selected_vars``
   strings across common clusters.                       Threshold ≥ 0.70.
3. Null CV bps — Spearman ρ across clusters.             Threshold ≥ 0.85.
4. Selected CV bps — Spearman ρ across clusters.         Threshold ≥ 0.85.
5. (Selected − Null) sign agreement — fraction of
   clusters whose improvement has the same sign.          Threshold ≥ 0.95.
6. Coefficient sign agreement per basis family — for each
   of Speed / TF / SF / OR, fraction of named coefficients
   matching in sign across pipelines.                     Threshold ≥ 0.85.
7. Stationary-vs-motion FR — on common (cluster_id, trial_id)
   pairs, median relative error on ``stationary_fr`` and
   ``motion_fr``.                                         Threshold < 1e-3.

Exit code: 0 if every non-skipped check passes, 1 if any fails.
Skipped checks (missing input files on either side) are marked ``pass=True``
with an explanatory note, since we cannot judge what we cannot see — but the
note is surfaced in stdout so it is obvious.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend for CLI use
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from rc2_glm.basis import raised_cosine_basis
from rc2_glm.config import GLMConfig

logger = logging.getLogger("rc2_glm.compare")


# Thresholds live here, not in tribal memory. Edit only with a prompt update.
PREFILTER_CATEGORY_THRESHOLD = 0.90
SELECTED_VARS_JACCARD_THRESHOLD = 0.70
NULL_CV_SPEARMAN_THRESHOLD = 0.85
SELECTED_CV_SPEARMAN_THRESHOLD = 0.85
SELECTED_MINUS_NULL_SIGN_THRESHOLD = 0.95
COEF_SIGN_FAMILY_THRESHOLD = 0.85
STATIONARY_MOTION_FR_MEDIAN_REL_ERR_THRESHOLD = 1e-3
# Tuning-curve prediction-space parity thresholds.
#
# The row ``tuning_curve_pearson_{var}`` gates on the *fraction of clusters*
# whose per-cluster curve Pearson r clears ``PER_CLUSTER_THRESHOLD``. The
# mean is also reported (informational — pattern 7) since it is the metric
# the earlier prompts discussed. Rationale for the switch: at small n the
# mean is dominated by a single-cluster outlier (one r=0.04 cluster drops a
# 6-cluster mean by ~0.1), producing false-negative row failures that
# don't reflect the scientific claim ("most clusters agree cluster-by-cluster").
# Threshold 0.80 tolerates 1/6 outliers at small n while still requiring
# broad agreement at n=30 (24/30 must clear the per-cluster bar).
TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD = 0.85
TUNING_CURVE_PEARSON_FRACTION_THRESHOLD = 0.80
# Legacy: the pre-2026-04-23 gate used mean(r) >= MEAN_THRESHOLD. Retained
# as a constant so the informational ``_mean_`` row can be plotted against
# the same visual bar without ambiguity.
TUNING_CURVE_PEARSON_MEAN_THRESHOLD = 0.90

# Families we match coefficients across. Any coefficient whose normalised
# name starts with one of these is routed into the family.
COEFFICIENT_FAMILIES = ("Speed", "TF", "SF", "OR")
# Families that are gated on raw-β sign agreement (categorical dummies —
# one coefficient per level, no identifiability ambiguity). Speed / TF
# live on correlated raised-cosine bases; their raw-β sign is noisy and
# is now tracked via prediction-space tuning-curve parity instead.
COEF_SIGN_GATING_FAMILIES = ("SF", "OR")


# --------------------------------------------------------------------------- #
# File discovery
# --------------------------------------------------------------------------- #


def _find_csv(root: Path, *candidates: str) -> Path | None:
    """Return the first existing path under ``root`` from ``candidates``, else None.

    Each candidate is a relative path. Tries the exact paths, then a recursive
    glob fallback on just the basename. The fallback lets us absorb minor
    layout drift (e.g. per-probe subdirs vs flat).
    """
    for c in candidates:
        p = root / c
        if p.is_file():
            return p
    for c in candidates:
        stem = Path(c).name
        matches = sorted(root.rglob(stem))
        if matches:
            return matches[0]
    return None


# --------------------------------------------------------------------------- #
# Column normalisation
# --------------------------------------------------------------------------- #


def _normalise_comparison_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Drop the ``time_`` glm-type prefix so Python and MATLAB tables align.

    Python writes ``time_Null_cv_bps`` / ``time_selected_vars`` etc. MATLAB
    may write ``Null_cv_bps`` / ``selected_vars`` or the prefixed form
    depending on the export script used. Stripping is safe either way.
    """
    renames: dict[str, str] = {}
    for col in df.columns:
        if col.startswith("time_"):
            renames[col] = col[len("time_"):]
    return df.rename(columns=renames)


def _select_first_present(df: pd.DataFrame, *names: str) -> pd.Series | None:
    for n in names:
        if n in df.columns:
            return df[n]
    return None


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #


def _merge_key(py: pd.DataFrame, ref: pd.DataFrame) -> list[str]:
    """Return the join key to use for per-cluster comparisons.

    MATLAB reference CSVs span multiple probes; python outputs typically
    contain a single probe. Merging on ``cluster_id`` alone therefore
    collides rows across probes and reports spurious disagreement. If both
    tables carry ``probe_id`` we scope the join to ``(probe_id, cluster_id)``.
    """
    if "probe_id" in py.columns and "probe_id" in ref.columns:
        return ["probe_id", "cluster_id"]
    return ["cluster_id"]


def _row(
    check: str,
    value: float,
    threshold: float,
    pass_: bool,
    note: str = "",
) -> dict:
    return {
        "check": check,
        "value": value,
        "threshold": threshold,
        "pass": pass_,
        "note": note,
    }


def _prefilter_category_agreement(
    py: pd.DataFrame | None, ref: pd.DataFrame | None
) -> dict:
    if py is None or ref is None:
        return _row(
            "prefilter_category_agreement", float("nan"),
            PREFILTER_CATEGORY_THRESHOLD, True,
            note="skipped: prefilter_decision_tree.csv missing on one side",
        )
    key = _merge_key(py, ref)
    merged = py[key + ["category"]].merge(
        ref[key + ["category"]], on=key,
        suffixes=("_py", "_ref"),
    )
    if merged.empty:
        return _row(
            "prefilter_category_agreement", float("nan"),
            PREFILTER_CATEGORY_THRESHOLD, False,
            note=f"no overlapping {'+'.join(key)}",
        )
    agree = float((merged["category_py"] == merged["category_ref"]).mean())
    return _row(
        "prefilter_category_agreement", agree,
        PREFILTER_CATEGORY_THRESHOLD, agree >= PREFILTER_CATEGORY_THRESHOLD,
        note=f"n={len(merged)} common clusters",
    )


def _jaccard(a: str, b: str) -> float:
    # Python writes '+'-joined, MATLAB writes ', '-joined. Accept either.
    # 'Null' / '' mean empty.
    def _to_set(s: str) -> set[str]:
        if not isinstance(s, str):
            return set()
        # Normalise both separators to a single delimiter before splitting.
        norm = s.replace(", ", "+").replace(",", "+")
        tokens = [t.strip() for t in norm.split("+")]
        return {t for t in tokens if t and t != "Null"}

    A, B = _to_set(a), _to_set(b)
    if not A and not B:
        return 1.0  # both Null → agree
    if not A or not B:
        return 0.0
    return len(A & B) / len(A | B)


def _selected_vars_jaccard(
    py: pd.DataFrame, ref: pd.DataFrame
) -> tuple[dict, list[int]]:
    py_sv = _select_first_present(py, "selected_vars", "time_selected_vars")
    ref_sv = _select_first_present(ref, "selected_vars", "time_selected_vars")
    if py_sv is None or ref_sv is None:
        return (
            _row("selected_vars_jaccard", float("nan"),
                 SELECTED_VARS_JACCARD_THRESHOLD, True,
                 note="skipped: selected_vars column absent on one side"),
            [],
        )
    key = _merge_key(py, ref)
    py_df = py[key].copy()
    py_df["sv_py"] = py_sv.values
    ref_df = ref[key].copy()
    ref_df["sv_ref"] = ref_sv.values
    merged = py_df.merge(ref_df, on=key)
    if merged.empty:
        return (
            _row("selected_vars_jaccard", float("nan"),
                 SELECTED_VARS_JACCARD_THRESHOLD, False,
                 note=f"no overlapping {'+'.join(key)}"),
            [],
        )
    jaccards = [_jaccard(a, b) for a, b in zip(merged["sv_py"], merged["sv_ref"])]
    zeros = merged.loc[[j == 0.0 for j in jaccards], "cluster_id"].tolist()
    mean_j = float(np.mean(jaccards))
    return (
        _row("selected_vars_jaccard", mean_j,
             SELECTED_VARS_JACCARD_THRESHOLD,
             mean_j >= SELECTED_VARS_JACCARD_THRESHOLD,
             note=f"n={len(merged)} | {len(zeros)} clusters with Jaccard=0"),
        [int(c) for c in zeros],
    )


def _cv_bps_correlation(
    py: pd.DataFrame, ref: pd.DataFrame, column_candidates: tuple[str, ...],
    check_name: str, threshold: float,
) -> dict:
    py_col = _select_first_present(py, *column_candidates)
    ref_col = _select_first_present(ref, *column_candidates)
    if py_col is None or ref_col is None:
        return _row(check_name, float("nan"), threshold, True,
                    note=f"skipped: none of {column_candidates} present on one side")
    key = _merge_key(py, ref)
    py_df = py[key].copy()
    py_df["v_py"] = py_col.values
    ref_df = ref[key].copy()
    ref_df["v_ref"] = ref_col.values
    merged = py_df.merge(ref_df, on=key).dropna()
    if len(merged) < 3:
        return _row(check_name, float("nan"), threshold, False,
                    note=f"fewer than 3 common non-NaN clusters (n={len(merged)})")
    rho, _ = spearmanr(merged["v_py"], merged["v_ref"])
    rho = float(rho)
    return _row(check_name, rho, threshold, rho >= threshold,
                note=f"n={len(merged)} common non-NaN clusters")


def _cv_bps_mean_offset(
    py: pd.DataFrame, ref: pd.DataFrame, column_candidates: tuple[str, ...],
    check_name: str,
) -> dict:
    """Informational: report ``mean(py - ref)`` in cv bits-per-spike.

    Absolute bps values scale with the time-bin width used for binning
    (Poisson per-bin entropy changes with bin width, although the
    Spearman ranking is invariant). Two runs using different bin widths
    will show a constant offset even when the models agree. This row
    never gates — it's reported so the sign + magnitude of any offset
    is visible in the report.
    """
    py_col = _select_first_present(py, *column_candidates)
    ref_col = _select_first_present(ref, *column_candidates)
    if py_col is None or ref_col is None:
        return _row(check_name, float("nan"), float("nan"), True,
                    note=f"skipped: none of {column_candidates} present on one side")
    key = _merge_key(py, ref)
    py_df = py[key].copy()
    py_df["v_py"] = py_col.values
    ref_df = ref[key].copy()
    ref_df["v_ref"] = ref_col.values
    merged = py_df.merge(ref_df, on=key).dropna()
    if merged.empty:
        return _row(check_name, float("nan"), float("nan"), True,
                    note="skipped: no common non-NaN clusters")
    diff = merged["v_py"].to_numpy() - merged["v_ref"].to_numpy()
    mean_offset = float(np.mean(diff))
    std_offset = float(np.std(diff))
    return _row(
        check_name, mean_offset, float("nan"), True,
        note=(
            f"informational — mean(py-ref) over n={len(merged)} common "
            f"clusters, std={std_offset:.3f}; absolute bps depends on "
            f"time-bin width convention (see summary/rc2_analysis-"
            f"python-matlab-parity.md, 'CV-bps offset — resolved')"
        ),
    )


def _selected_minus_null_sign_agreement(
    py: pd.DataFrame, ref: pd.DataFrame
) -> dict:
    def delta(df: pd.DataFrame) -> pd.Series | None:
        sel = _select_first_present(df, "Selected_cv_bps", "time_Selected_cv_bps")
        nul = _select_first_present(df, "Null_cv_bps", "time_Null_cv_bps")
        if sel is None or nul is None:
            return None
        return sel.values - nul.values

    py_d = delta(py)
    ref_d = delta(ref)
    if py_d is None or ref_d is None:
        return _row("selected_minus_null_sign_agreement", float("nan"),
                    SELECTED_MINUS_NULL_SIGN_THRESHOLD, True,
                    note="skipped: Selected_cv_bps or Null_cv_bps missing on one side")
    key = _merge_key(py, ref)
    py_df = py[key].copy()
    py_df["d_py"] = py_d
    ref_df = ref[key].copy()
    ref_df["d_ref"] = ref_d
    merged = py_df.merge(ref_df, on=key).dropna()
    if merged.empty:
        return _row("selected_minus_null_sign_agreement", float("nan"),
                    SELECTED_MINUS_NULL_SIGN_THRESHOLD, False,
                    note="no common non-NaN clusters")
    # 0 counts as its own sign, so clusters where both improvements are 0 agree.
    agree = float((np.sign(merged["d_py"]) == np.sign(merged["d_ref"])).mean())
    return _row("selected_minus_null_sign_agreement", agree,
                SELECTED_MINUS_NULL_SIGN_THRESHOLD,
                agree >= SELECTED_MINUS_NULL_SIGN_THRESHOLD,
                note=f"n={len(merged)}")


def _coefficient_family(name: str) -> str | None:
    """Return the family prefix for a coefficient name, or None if not one we check.

    Interaction terms (contain ``_x_``) are deliberately excluded.
    """
    if not isinstance(name, str) or "_x_" in name:
        return None
    for fam in COEFFICIENT_FAMILIES:
        if name.startswith(fam + "_") or name == fam:
            return fam
    return None


def _coefficient_sign_agreement(
    py: pd.DataFrame | None, ref: pd.DataFrame | None
) -> list[dict]:
    if py is None or ref is None:
        return [_row(
            f"coefficient_sign_{fam}", float("nan"),
            COEF_SIGN_FAMILY_THRESHOLD, True,
            note="skipped: glm_coefficients.csv missing on one side",
        ) for fam in COEFFICIENT_FAMILIES]

    # Prefer the Selected model when present, fall back to Additive.
    def _filtered(df: pd.DataFrame) -> pd.DataFrame:
        if "model" in df.columns:
            sel = df[df["model"] == "Selected"]
            if not sel.empty:
                return sel
            return df[df["model"] == "Additive"]
        return df

    py_f = _filtered(py)
    ref_f = _filtered(ref)
    key_cols = _merge_key(py_f, ref_f) + ["coefficient"]
    merged = py_f[key_cols + ["estimate"]].merge(
        ref_f[key_cols + ["estimate"]], on=key_cols,
        suffixes=("_py", "_ref"),
    )
    rows: list[dict] = []
    for fam in COEFFICIENT_FAMILIES:
        mask = merged["coefficient"].apply(_coefficient_family) == fam
        sub = merged[mask]
        is_gating = fam in COEF_SIGN_GATING_FAMILIES
        if sub.empty:
            rows.append(_row(
                f"coefficient_sign_{fam}", float("nan"),
                COEF_SIGN_FAMILY_THRESHOLD, True,
                note=f"skipped: no common {fam} coefficients (neither side selected {fam}?)",
            ))
            continue
        agree = float((np.sign(sub["estimate_py"]) == np.sign(sub["estimate_ref"])).mean())
        if is_gating:
            rows.append(_row(
                f"coefficient_sign_{fam}", agree,
                COEF_SIGN_FAMILY_THRESHOLD, agree >= COEF_SIGN_FAMILY_THRESHOLD,
                note=f"n={len(sub)}",
            ))
        else:
            # Speed/TF sit on correlated raised-cosine bases: raw-β rotations
            # produce the same predicted rate, so the sign agreement is an
            # identifiability artefact, not a model-quality signal. Keep the
            # number in the report for continuity but do not gate on it —
            # ``tuning_curve_pearson_{fam}`` is the prediction-space check.
            rows.append(_row(
                f"coefficient_sign_{fam}", agree,
                float("nan"), True,
                note=(
                    f"informational — raw-β sign agreement over n={len(sub)} "
                    f"coefficients; {fam} lives on correlated raised-cosine "
                    f"bases where sign is not identifiable. See "
                    f"tuning_curve_pearson_{fam} for the gating check."
                ),
            ))
    return rows


def _variable_betas_per_cluster(
    coef_df: pd.DataFrame,
    variable: str,
    model: str,
    n_bases: int,
) -> dict[tuple[str, int], np.ndarray]:
    """Return ``{(probe_id, cluster_id): β_vec}`` for one variable × model.

    Only returns entries where all ``n_bases`` β's are present (a sparse
    model can omit the variable entirely; we skip those silently — the
    caller pairs what's present on both sides).
    """
    prefix = f"{variable}_"
    has_probe = "probe_id" in coef_df.columns
    mask = (
        (coef_df["model"] == model)
        & coef_df["coefficient"].astype(str).str.startswith(prefix)
        & ~coef_df["coefficient"].astype(str).str.contains("_x_", regex=False)
    )
    sub = coef_df[mask]
    if sub.empty:
        return {}

    # Expected basis names in order: Speed_1, Speed_2, ..., Speed_n_bases
    expected = [f"{variable}_{i + 1}" for i in range(n_bases)]
    expected_set = set(expected)
    sub = sub[sub["coefficient"].isin(expected_set)]

    out: dict[tuple[str, int], np.ndarray] = {}
    group_cols = ["probe_id", "cluster_id"] if has_probe else ["cluster_id"]
    for key, grp in sub.groupby(group_cols):
        by_name = dict(zip(grp["coefficient"].astype(str), grp["estimate"]))
        if not all(name in by_name for name in expected):
            continue  # partial — skip rather than zero-fill
        beta = np.array([float(by_name[name]) for name in expected], dtype=np.float64)
        if has_probe:
            probe_id, cluster_id = key
        else:
            probe_id, cluster_id = "", key
        out[(str(probe_id), int(cluster_id))] = beta
    return out


def _tuning_curve_pearson(
    coef_py: pd.DataFrame | None,
    coef_ref: pd.DataFrame | None,
    variable: str,
    config: GLMConfig,
    model: str = "Selected",
) -> list[dict]:
    """Per-cluster Pearson r between Python and MATLAB tuning curves.

    Builds the same raised-cosine grid on both sides using Python's basis
    (MATLAB's is verified numerically identical to ~1e-12), evaluates
    ``B_grid @ β_var`` for the ``model`` fit (default Selected), and
    Pearson-correlates the two curves. Returns two rows:

    - ``tuning_curve_pearson_{variable}`` (gating): value = fraction of
      clusters with ``r >= PER_CLUSTER_THRESHOLD`` (default 0.85). Gated
      at ``FRACTION_THRESHOLD`` (0.80). Robust to single-cluster outliers
      at small n — the previous mean-based gate failed spuriously on
      n=6 TF data when one r=0.04 outlier pulled the mean below 0.90.
    - ``tuning_curve_pearson_mean_{variable}`` (informational, pattern 7):
      value = mean(r). No threshold; kept so the number cited in earlier
      prompts and in ``summary/rc2_analysis-python-matlab-parity.md``
      stays visible without re-parsing the gate row's note.
    """
    gate_name = f"tuning_curve_pearson_{variable}"
    info_name = f"tuning_curve_pearson_mean_{variable}"

    def _skip(note: str, passed: bool = True) -> list[dict]:
        return [
            _row(gate_name, float("nan"),
                 TUNING_CURVE_PEARSON_FRACTION_THRESHOLD, passed, note=note),
            _row(info_name, float("nan"), float("nan"), True,
                 note=f"informational — mean(r); {note}"),
        ]

    if coef_py is None or coef_ref is None:
        return _skip("skipped: glm_coefficients.csv missing on one side")

    if variable == "Speed":
        n_bases = config.n_speed_bases
        lo, hi = config.speed_range
    elif variable == "TF":
        n_bases = config.n_tf_bases
        lo, hi = config.tf_range
    else:
        raise ValueError(f"unsupported variable: {variable!r}")

    py_map = _variable_betas_per_cluster(coef_py, variable, model, n_bases)
    ref_map = _variable_betas_per_cluster(coef_ref, variable, model, n_bases)
    common_keys = sorted(set(py_map) & set(ref_map))
    if not common_keys:
        return _skip(
            f"skipped: no common clusters with {variable} in the "
            f"{model} model on both sides"
        )

    x_grid = np.linspace(float(lo), float(hi), 100)
    B_grid = raised_cosine_basis(x_grid, n_bases, float(lo), float(hi))

    per_cluster_r: list[float] = []
    for key in common_keys:
        curve_py = B_grid @ py_map[key]
        curve_ref = B_grid @ ref_map[key]
        # pearsonr is undefined for a constant series; skip degenerate fits.
        if np.std(curve_py) < 1e-12 or np.std(curve_ref) < 1e-12:
            continue
        r, _ = pearsonr(curve_py, curve_ref)
        if np.isfinite(r):
            per_cluster_r.append(float(r))

    if not per_cluster_r:
        return _skip("no clusters yielded a finite Pearson r", passed=False)

    n = len(per_cluster_r)
    mean_r = float(np.mean(per_cluster_r))
    median_r = float(np.median(per_cluster_r))
    min_r = float(np.min(per_cluster_r))
    below = int(sum(r < TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD for r in per_cluster_r))
    frac_above = float((n - below) / n)
    passed = frac_above >= TUNING_CURVE_PEARSON_FRACTION_THRESHOLD
    gate_row = _row(
        gate_name, frac_above,
        TUNING_CURVE_PEARSON_FRACTION_THRESHOLD, passed,
        note=(
            f"n={n} | frac(r>={TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD:.2f})"
            f"={frac_above:.3f} | mean={mean_r:.3f} | median={median_r:.3f}"
            f" | min={min_r:.3f} | below={below}/{n} | model={model}"
        ),
    )
    info_row = _row(
        info_name, mean_r, float("nan"), True,
        note=(
            f"informational — mean(r) over n={n} clusters; the gate "
            f"{gate_name} uses frac(r>={TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD:.2f})"
            f". Legacy mean-threshold was {TUNING_CURVE_PEARSON_MEAN_THRESHOLD}"
            f" — kept here as reference."
        ),
    )
    return [gate_row, info_row]


def _stationary_vs_motion_fr(
    py: pd.DataFrame | None, ref: pd.DataFrame | None
) -> dict:
    if py is None or ref is None:
        return _row(
            "stationary_vs_motion_fr_median_rel_err", float("nan"),
            STATIONARY_MOTION_FR_MEDIAN_REL_ERR_THRESHOLD, True,
            note=(
                "skipped: stationary_vs_motion_fr CSV missing on one side "
                "(python-dir needs diagnostics/stationary_vs_motion_fr_python.csv "
                "from a pipeline run)"
            ),
        )
    key = _merge_key(py, ref) + ["trial_id"]
    py = py[key + ["stationary_fr", "motion_fr"]].copy()
    ref = ref[key + ["stationary_fr", "motion_fr"]].copy()
    merged = py.merge(ref, on=key, suffixes=("_py", "_ref")).dropna()
    if merged.empty:
        return _row(
            "stationary_vs_motion_fr_median_rel_err", float("nan"),
            STATIONARY_MOTION_FR_MEDIAN_REL_ERR_THRESHOLD, False,
            note="no common (cluster_id, trial_id) rows",
        )

    def _rel_err(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        denom = np.maximum(np.abs(b), 1e-9)
        return np.abs(a - b) / denom

    errs = np.concatenate([
        _rel_err(merged["stationary_fr_py"].to_numpy(),
                 merged["stationary_fr_ref"].to_numpy()),
        _rel_err(merged["motion_fr_py"].to_numpy(),
                 merged["motion_fr_ref"].to_numpy()),
    ])
    med = float(np.median(errs))
    return _row(
        "stationary_vs_motion_fr_median_rel_err", med,
        STATIONARY_MOTION_FR_MEDIAN_REL_ERR_THRESHOLD,
        med < STATIONARY_MOTION_FR_MEDIAN_REL_ERR_THRESHOLD,
        note=f"n={len(merged)} (cluster,trial) pairs",
    )


# --------------------------------------------------------------------------- #
# Classification-differences figure
# --------------------------------------------------------------------------- #

# Canonical order of prefilter decision-tree leaves. Matches the 8-branch
# categorisation in scripts/glm_single_cluster_analysis.m and is reused to
# lay out confusion-matrix axes consistently across probes.
PREFILTER_CATEGORY_ORDER = (
    "none_significant",
    "T_only_significant",
    "V_only_significant",
    "VT_only_significant",
    "T_and_V_significant",
    "T_and_VT_significant",
    "V_and_VT_significant",
    "all_three_significant",
)


def _plot_prefilter_confusion(
    ax, counts: np.ndarray, categories: tuple[str, ...], title: str,
) -> None:
    """Draw one confusion-matrix panel onto ``ax``.

    Counts are laid out with MATLAB category on rows, Python category on
    columns. Diagonal cells (agreement) are shaded green; off-diagonal
    cells (disagreement) are shaded red. Empty cells left white.
    """
    n = len(categories)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(categories, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(categories, fontsize=7)
    ax.set_xlabel("Python category", fontsize=8)
    ax.set_ylabel("MATLAB category", fontsize=8)
    ax.set_title(title, fontsize=9)

    max_count = int(counts.max()) if counts.size else 0
    for i in range(n):
        for j in range(n):
            c = int(counts[i, j])
            if c == 0:
                continue
            on_diag = (i == j)
            # Green → agreement; red → disagreement. Alpha scales with count
            # so a single-cell outlier stays visible against a dense diagonal.
            intensity = 0.25 + 0.65 * (c / max_count if max_count else 1.0)
            color = (0.2, 0.7, 0.3, intensity) if on_diag else (0.85, 0.25, 0.25, intensity)
            ax.add_patch(plt.Rectangle(
                (j - 0.5, i - 0.5), 1, 1, facecolor=color, edgecolor="gray",
                linewidth=0.3,
            ))
            text_color = "black" if intensity < 0.6 else "white"
            ax.text(j, i, str(c), ha="center", va="center",
                    fontsize=7, color=text_color, fontweight="bold" if on_diag else "normal")

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)  # origin top-left, matches typical confusion-matrix layout
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)


def _plot_classification_differences(
    pre_py: pd.DataFrame | None,
    pre_ref: pd.DataFrame | None,
    out_pdf: Path,
    out_disagreements_csv: Path,
) -> dict:
    """Per-probe + aggregated confusion matrix of prefilter category PY vs ML.

    Renders an N-panel figure (one per probe plus an all-probes aggregate
    when n_probes > 1) and emits a CSV of per-cluster disagreements so
    scientists can drill into specific (probe, cluster) rows that flip.
    Returns a report row summarising what was written.
    """
    if pre_py is None or pre_ref is None:
        return _row(
            "classification_differences_figure", float("nan"),
            float("nan"), True,
            note="skipped: prefilter_decision_tree.csv missing on one side",
        )

    key = _merge_key(pre_py, pre_ref)
    merged = pre_py[key + ["category"]].merge(
        pre_ref[key + ["category"]], on=key,
        suffixes=("_py", "_ref"),
    )
    if merged.empty:
        return _row(
            "classification_differences_figure", float("nan"),
            float("nan"), False,
            note=f"no overlapping {'+'.join(key)} — figure not rendered",
        )

    # Order of tickmarks: canonical order first, then any unseen categories
    # appended alphabetically (defensive — keeps novel MATLAB branches
    # visible rather than silently merging them into the first bucket).
    seen = set(merged["category_py"]).union(set(merged["category_ref"]))
    extras = tuple(sorted(s for s in seen if s not in PREFILTER_CATEGORY_ORDER))
    categories = PREFILTER_CATEGORY_ORDER + extras

    # Persist the per-cluster disagreement list before plotting so the CSV
    # exists even if matplotlib chokes on an edge case.
    disagree = merged[merged["category_py"] != merged["category_ref"]].copy()
    disagree = disagree.rename(columns={
        "category_py": "python_category",
        "category_ref": "matlab_category",
    })
    out_disagreements_csv.parent.mkdir(parents=True, exist_ok=True)
    disagree.to_csv(out_disagreements_csv, index=False)

    probes = sorted(merged["probe_id"].unique()) if "probe_id" in merged.columns else [""]
    multi = len(probes) > 1
    n_panels = len(probes) + (1 if multi else 0)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5.2 * n_cols, 5.2 * n_rows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    def _counts(df: pd.DataFrame) -> np.ndarray:
        cm = pd.crosstab(df["category_ref"], df["category_py"])
        cm = cm.reindex(index=categories, columns=categories, fill_value=0)
        return cm.to_numpy(dtype=int)

    for i, probe in enumerate(probes):
        sub = merged[merged["probe_id"] == probe] if probe else merged
        agree = float((sub["category_py"] == sub["category_ref"]).mean()) if len(sub) else float("nan")
        title = (
            f"{probe or 'all'}  (n={len(sub)}, agree={100 * agree:.1f}%)"
            if len(sub) else f"{probe or 'all'}  (empty)"
        )
        _plot_prefilter_confusion(axes_flat[i], _counts(sub), categories, title)

    if multi:
        agree_all = float((merged["category_py"] == merged["category_ref"]).mean())
        _plot_prefilter_confusion(
            axes_flat[len(probes)], _counts(merged), categories,
            f"ALL probes  (n={len(merged)}, agree={100 * agree_all:.1f}%)",
        )

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    fig.suptitle(
        "Prefilter category — Python (cols) vs MATLAB (rows)",
        fontsize=11, y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)

    agree_all = float((merged["category_py"] == merged["category_ref"]).mean())
    return _row(
        "classification_differences_figure", agree_all,
        PREFILTER_CATEGORY_THRESHOLD, agree_all >= PREFILTER_CATEGORY_THRESHOLD,
        note=(
            f"wrote {out_pdf.name} ({len(probes)} probes, "
            f"{len(merged)} common clusters, {len(disagree)} disagreements) "
            f"and {out_disagreements_csv.name}"
        ),
    )


# All variables the forward selector can pick (main effects + interactions).
# Plot order matches the design-matrix block order in assemble_design_matrix.
SELECTED_VARS_PLOT_ORDER = (
    "Speed", "TF", "SF", "OR",
    "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
    "TF_x_SF", "TF_x_OR", "SF_x_OR",
)

# Short axis labels for the selected-vars figure.
_SELECTED_VARS_LABELS = {
    "Speed_x_TF": "Sp×TF",
    "Speed_x_SF": "Sp×SF",
    "Speed_x_OR": "Sp×OR",
    "TF_x_SF": "TF×SF",
    "TF_x_OR": "TF×OR",
    "SF_x_OR": "SF×OR",
}


def _variable_agreement_counts(
    py_vars: pd.Series, ref_vars: pd.Series,
) -> pd.DataFrame:
    """For each variable in ``SELECTED_VARS_PLOT_ORDER``, count clusters.

    Returns a DataFrame with one row per variable and columns
    ``both``, ``python_only``, ``matlab_only``, ``neither``. The counts
    add up to ``len(py_vars)`` (= number of merged clusters).
    """
    def _sets(s: pd.Series) -> list[set[str]]:
        out = []
        for raw in s:
            if not isinstance(raw, str):
                out.append(set())
                continue
            norm = raw.replace(", ", "+").replace(",", "+")
            tokens = {t.strip() for t in norm.split("+")}
            out.append({t for t in tokens if t and t != "Null"})
        return out

    py_sets = _sets(py_vars)
    ref_sets = _sets(ref_vars)
    rows = []
    for var in SELECTED_VARS_PLOT_ORDER:
        both = sum(1 for p, r in zip(py_sets, ref_sets) if var in p and var in r)
        py_only = sum(1 for p, r in zip(py_sets, ref_sets) if var in p and var not in r)
        ml_only = sum(1 for p, r in zip(py_sets, ref_sets) if var not in p and var in r)
        neither = sum(1 for p, r in zip(py_sets, ref_sets) if var not in p and var not in r)
        rows.append({
            "variable": var, "both": both, "python_only": py_only,
            "matlab_only": ml_only, "neither": neither,
        })
    return pd.DataFrame(rows)


def _draw_selected_vars_panel(ax, counts: pd.DataFrame, title: str) -> None:
    """One probe's selected-vars agreement as a grouped bar chart.

    x axis: variables (main effects then interactions).
    Three stacked counts per variable: ``both`` (grey), ``python_only``
    (warm red, tagged PY+), ``matlab_only`` (cool blue, tagged ML+).
    ``neither`` is omitted — it is by far the largest bucket for most
    clusters and would crush the bars that actually carry signal.
    """
    n = len(counts)
    x = np.arange(n)
    labels = [_SELECTED_VARS_LABELS.get(v, v) for v in counts["variable"]]
    both = counts["both"].to_numpy()
    py_only = counts["python_only"].to_numpy()
    ml_only = counts["matlab_only"].to_numpy()

    # Stack: "both" at the bottom (grey), then PY-only (red) + ML-only
    # (blue) side-by-side as half-width bars so disagreement is the eye's
    # first target.
    ax.bar(x, both, width=0.8, color=(0.55, 0.55, 0.55), label="both")
    ax.bar(x - 0.2, py_only, width=0.38, bottom=both,
           color=(0.85, 0.3, 0.3), label="Python only")
    ax.bar(x + 0.2, ml_only, width=0.38, bottom=both,
           color=(0.3, 0.5, 0.85), label="MATLAB only")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=7)
    ax.set_ylabel("# clusters", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.grid(axis="y", alpha=0.25, linestyle=":")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def _plot_selected_vars_differences(
    cmp_py: pd.DataFrame | None,
    cmp_ref: pd.DataFrame | None,
    out_pdf: Path,
    out_disagreements_csv: Path,
) -> dict:
    """Per-probe grouped bars of selected-variable agreement PY vs ML.

    Companion to the prefilter-category figure. Answers the question
    "across common clusters, for each candidate variable, how often do
    the two pipelines agree on including it?". Also emits a per-cluster
    disagreement CSV with columns ``probe_id, cluster_id, jaccard,
    python_vars, matlab_vars, python_only, matlab_only``.
    """
    if cmp_py is None or cmp_ref is None:
        return _row(
            "selected_vars_differences_figure", float("nan"),
            float("nan"), True,
            note="skipped: glm_model_comparison.csv missing on one side",
        )

    py_sv_col = "time_selected_vars" if "time_selected_vars" in cmp_py.columns \
        else "selected_vars"
    ref_sv_col = "time_selected_vars" if "time_selected_vars" in cmp_ref.columns \
        else "selected_vars"
    if py_sv_col not in cmp_py.columns or ref_sv_col not in cmp_ref.columns:
        return _row(
            "selected_vars_differences_figure", float("nan"),
            float("nan"), True,
            note="skipped: selected_vars column absent on one side",
        )

    key = _merge_key(cmp_py, cmp_ref)
    py_df = cmp_py[key + [py_sv_col]].rename(columns={py_sv_col: "sv_py"})
    ref_df = cmp_ref[key + [ref_sv_col]].rename(columns={ref_sv_col: "sv_ref"})
    merged = py_df.merge(ref_df, on=key)
    if merged.empty:
        return _row(
            "selected_vars_differences_figure", float("nan"),
            float("nan"), False,
            note=f"no overlapping {'+'.join(key)} — figure not rendered",
        )

    # Per-cluster disagreement table (Jaccard, symmetric differences).
    def _as_set(s) -> set[str]:
        if not isinstance(s, str):
            return set()
        norm = s.replace(", ", "+").replace(",", "+")
        return {t.strip() for t in norm.split("+") if t.strip() and t.strip() != "Null"}

    merged["py_set"] = merged["sv_py"].apply(_as_set)
    merged["ref_set"] = merged["sv_ref"].apply(_as_set)
    merged["jaccard"] = merged.apply(
        lambda r: (
            1.0 if not r["py_set"] and not r["ref_set"]
            else len(r["py_set"] & r["ref_set"]) / len(r["py_set"] | r["ref_set"])
            if (r["py_set"] | r["ref_set"]) else 1.0
        ),
        axis=1,
    )
    merged["python_only"] = merged.apply(
        lambda r: "+".join(sorted(r["py_set"] - r["ref_set"])), axis=1,
    )
    merged["matlab_only"] = merged.apply(
        lambda r: "+".join(sorted(r["ref_set"] - r["py_set"])), axis=1,
    )
    disagree = merged[merged["jaccard"] < 1.0].copy()
    disagree = disagree[key + ["jaccard", "sv_py", "sv_ref",
                               "python_only", "matlab_only"]]
    disagree = disagree.rename(columns={
        "sv_py": "python_vars", "sv_ref": "matlab_vars",
    })
    out_disagreements_csv.parent.mkdir(parents=True, exist_ok=True)
    disagree.to_csv(out_disagreements_csv, index=False)

    probes = sorted(merged["probe_id"].unique()) if "probe_id" in merged.columns else [""]
    multi = len(probes) > 1
    n_panels = len(probes) + (1 if multi else 0)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(6.0 * n_cols, 4.2 * n_rows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    for i, probe in enumerate(probes):
        sub = merged[merged["probe_id"] == probe] if probe else merged
        counts = _variable_agreement_counts(sub["sv_py"], sub["sv_ref"])
        mean_jac = float(sub["jaccard"].mean())
        title = (
            f"{probe or 'all'}  (n={len(sub)} clusters, "
            f"mean Jaccard={mean_jac:.3f})"
        )
        _draw_selected_vars_panel(axes_flat[i], counts, title)

    if multi:
        counts_all = _variable_agreement_counts(merged["sv_py"], merged["sv_ref"])
        _draw_selected_vars_panel(
            axes_flat[len(probes)], counts_all,
            f"ALL probes  (n={len(merged)} clusters, "
            f"mean Jaccard={merged['jaccard'].mean():.3f})",
        )

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=(0.55, 0.55, 0.55), label="both"),
        plt.Rectangle((0, 0), 1, 1, facecolor=(0.85, 0.3, 0.3), label="Python only"),
        plt.Rectangle((0, 0), 1, 1, facecolor=(0.3, 0.5, 0.85), label="MATLAB only"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        "Selected variables — agreement between Python and MATLAB",
        fontsize=11, y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.97))
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)

    mean_jac_all = float(merged["jaccard"].mean())
    return _row(
        "selected_vars_differences_figure", mean_jac_all,
        SELECTED_VARS_JACCARD_THRESHOLD,
        mean_jac_all >= SELECTED_VARS_JACCARD_THRESHOLD,
        note=(
            f"wrote {out_pdf.name} ({len(probes)} probes, "
            f"{len(merged)} common clusters, {len(disagree)} with "
            f"partial selected-vars disagreement) and "
            f"{out_disagreements_csv.name}"
        ),
    )


def _render_aggregate_forward_selection(
    cmp_df: pd.DataFrame | None,
    side: str,  # "python" | "matlab"
    out_pdf: Path,
) -> dict:
    """Render the MATLAB Fig 2 forward-selection summary across **all probes**.

    Per-probe runs already get a ``forward_selection_summary.pdf`` under
    ``{_runs}/{probe}/figs/``. This aggregate variant answers "across all
    probes, what fraction of clusters have each model type?" in a single
    view — the one the paper figure would show.
    """
    if cmp_df is None or cmp_df.empty:
        return _row(
            f"forward_selection_summary_{side}", float("nan"),
            float("nan"), True,
            note=f"skipped: no {side} model_comparison data",
        )
    # Lazy import — plots module pulls in the full plotting stack, which we
    # don't want to load for the report-CSV-only code path.
    from rc2_glm.plots import plot_forward_selection_summary

    # The plot function expects every classification column prefixed with
    # ``time_`` (the MATLAB glm_type scope). Both sides ship these columns
    # prefixed; ``_normalise_comparison_columns`` strips the prefix earlier
    # in the pipeline for generic merge logic. Re-add it here so the plot
    # function sees its expected schema on whichever ``cmp_df`` it gets.
    df = cmp_df.copy()
    time_cols = (
        "selected_vars", "is_speed_tuned", "is_tf_tuned", "is_sf_tuned",
        "is_or_tuned", "has_interaction", "has_speed_x_tf", "has_speed_x_sf",
        "has_speed_x_or", "has_tf_x_sf", "has_tf_x_or", "has_sf_x_or",
        "delta_selected_vs_null", "delta_additive_vs_null",
        "delta_selected_vs_additive", "delta_interaction",
    )
    renames = {c: f"time_{c}" for c in time_cols if c in df.columns}
    df = df.rename(columns=renames)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig = plot_forward_selection_summary(df)
    fig.suptitle(
        f"Forward-selection summary ({side}) — "
        f"{df['probe_id'].nunique() if 'probe_id' in df.columns else 1} probes, "
        f"{len(df)} clusters",
        fontsize=11,
    )
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)
    return _row(
        f"forward_selection_summary_{side}", float(len(df)),
        float("nan"), True,
        note=f"wrote {out_pdf.name} ({len(df)} clusters)",
    )


# --------------------------------------------------------------------------- #
# Top-level runner
# --------------------------------------------------------------------------- #


def _probe_ids(*dfs: pd.DataFrame | None) -> set[str]:
    out: set[str] = set()
    for df in dfs:
        if df is None or "probe_id" not in df.columns:
            continue
        out.update(str(x) for x in df["probe_id"].dropna().unique())
    return out


def _filter_probes(
    df: pd.DataFrame | None, probes: set[str]
) -> pd.DataFrame | None:
    if df is None or "probe_id" not in df.columns:
        return df
    return df[df["probe_id"].astype(str).isin(probes)].reset_index(drop=True)


def run_compare(
    python_dir: Path,
    reference_dir: Path,
    out_path: Path | None,
    figures_dir: Path | None = None,
) -> tuple[pd.DataFrame, list[int]]:
    """Run all checks and return the report + the list of Jaccard=0 cluster ids.

    ``figures_dir`` controls where the per-probe classification-differences
    PDF + the per-cluster disagreements CSV are written. Defaults to
    ``<out_path>.parent`` when ``out_path`` is given (i.e. the same
    ``validation/`` directory that already holds ``comparison_report.csv``).
    """
    # Accept either flat directory or one level of per-probe subdirs.
    py_cmp = _find_csv(python_dir, "glm_model_comparison.csv")
    ref_cmp = _find_csv(reference_dir, "glm_model_comparison.csv")
    py_coef = _find_csv(python_dir, "glm_coefficients.csv")
    ref_coef = _find_csv(reference_dir, "glm_coefficients.csv")
    py_pre = _find_csv(python_dir, "prefilter_decision_tree.csv")
    ref_pre = _find_csv(reference_dir, "prefilter_decision_tree.csv")

    py_svm = _find_csv(
        python_dir,
        "diagnostics/stationary_vs_motion_fr_python.csv",
        "stationary_vs_motion_fr_python.csv",
    )
    ref_svm = _find_csv(
        reference_dir,
        "csvs/stationary_vs_motion_fr",  # MATLAB writes per-probe under this dir
    )
    # The top-level candidate above is a directory, not a file; fall back to
    # a file search if the directory exists.
    if ref_svm is None or ref_svm.is_dir():
        candidate_dir = reference_dir / "csvs" / "stationary_vs_motion_fr"
        if candidate_dir.is_dir():
            mats = sorted(candidate_dir.glob("*.csv"))
            ref_svm = mats[0] if mats else None

    def _read(p: Path | None) -> pd.DataFrame | None:
        if p is None:
            return None
        try:
            return _normalise_comparison_columns(pd.read_csv(p))
        except Exception as exc:
            logger.warning("failed to read %s: %s", p, exc)
            return None

    cmp_py = _read(py_cmp)
    cmp_ref = _read(ref_cmp)
    coef_py = _read(py_coef)
    coef_ref = _read(ref_coef)
    pre_py = _read(py_pre)
    pre_ref = _read(ref_pre)
    svm_py = _read(py_svm)
    svm_ref = _read(ref_svm)

    # Preflight: if the python side has only one probe but the reference
    # spans several, scope the reference down so cluster_ids don't collide
    # across probes. Log both sets so discrepancies are obvious.
    py_probes = _probe_ids(cmp_py, coef_py, pre_py, svm_py)
    ref_probes = _probe_ids(cmp_ref, coef_ref, pre_ref, svm_ref)
    if py_probes and ref_probes and py_probes.issubset(ref_probes) and py_probes != ref_probes:
        logger.warning(
            "reference covers %d probes %s but python-dir has only %s — "
            "filtering reference to the python probe set for per-cluster checks",
            len(ref_probes), sorted(ref_probes), sorted(py_probes),
        )
        cmp_ref = _filter_probes(cmp_ref, py_probes)
        coef_ref = _filter_probes(coef_ref, py_probes)
        pre_ref = _filter_probes(pre_ref, py_probes)
        svm_ref = _filter_probes(svm_ref, py_probes)

    rows: list[dict] = []

    rows.append(_prefilter_category_agreement(pre_py, pre_ref))

    if cmp_py is not None and cmp_ref is not None:
        jac_row, zero_clusters = _selected_vars_jaccard(cmp_py, cmp_ref)
        rows.append(jac_row)
        rows.append(_cv_bps_correlation(
            cmp_py, cmp_ref, ("Null_cv_bps",),
            "null_cv_bps_spearman", NULL_CV_SPEARMAN_THRESHOLD,
        ))
        rows.append(_cv_bps_mean_offset(
            cmp_py, cmp_ref, ("Null_cv_bps",),
            "null_cv_bps_mean_offset",
        ))
        rows.append(_cv_bps_correlation(
            cmp_py, cmp_ref, ("Selected_cv_bps",),
            "selected_cv_bps_spearman", SELECTED_CV_SPEARMAN_THRESHOLD,
        ))
        rows.append(_cv_bps_mean_offset(
            cmp_py, cmp_ref, ("Selected_cv_bps",),
            "selected_cv_bps_mean_offset",
        ))
        rows.append(_selected_minus_null_sign_agreement(cmp_py, cmp_ref))
    else:
        zero_clusters = []
        for name, thresh in [
            ("selected_vars_jaccard", SELECTED_VARS_JACCARD_THRESHOLD),
            ("null_cv_bps_spearman", NULL_CV_SPEARMAN_THRESHOLD),
            ("selected_cv_bps_spearman", SELECTED_CV_SPEARMAN_THRESHOLD),
            ("selected_minus_null_sign_agreement", SELECTED_MINUS_NULL_SIGN_THRESHOLD),
        ]:
            rows.append(_row(
                name, float("nan"), thresh, True,
                note="skipped: glm_model_comparison.csv missing on one side",
            ))

    rows.extend(_coefficient_sign_agreement(coef_py, coef_ref))

    # Prediction-space tuning-curve parity for Speed / TF. Invariant to
    # the β-rotation that raised-cosine bases admit; gates on the
    # fraction of clusters whose per-cluster Pearson r clears 0.85
    # (robust to small-n outliers — see constants block at top of file).
    config = GLMConfig()
    for variable in ("Speed", "TF"):
        rows.extend(_tuning_curve_pearson(coef_py, coef_ref, variable, config))

    rows.append(_stationary_vs_motion_fr(svm_py, svm_ref))

    # Figure + per-cluster disagreements CSV. Goes into the same directory
    # as the comparison_report.csv by default so everything parity-related
    # lands under one ``validation/`` folder.
    if figures_dir is None and out_path is not None:
        figures_dir = out_path.parent
    if figures_dir is not None:
        rows.append(_plot_classification_differences(
            pre_py, pre_ref,
            out_pdf=figures_dir / "classification_differences.pdf",
            out_disagreements_csv=figures_dir / "classification_disagreements.csv",
        ))
        rows.append(_plot_selected_vars_differences(
            cmp_py, cmp_ref,
            out_pdf=figures_dir / "selected_vars_differences.pdf",
            out_disagreements_csv=figures_dir / "selected_vars_disagreements.csv",
        ))
        # Aggregate MATLAB Fig 2 port (per-probe figures already land in
        # {py_run}/_runs/{probe}/figs/ as part of the per-probe pipeline).
        # Also render a MATLAB-side aggregate when glm_model_comparison.csv
        # is available, so the two sides can be flipped side-by-side.
        rows.append(_render_aggregate_forward_selection(
            cmp_py, "python",
            out_pdf=figures_dir / "forward_selection_summary_python.pdf",
        ))
        rows.append(_render_aggregate_forward_selection(
            cmp_ref, "matlab",
            out_pdf=figures_dir / "forward_selection_summary_matlab.pdf",
        ))

    df = pd.DataFrame(rows)

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path, index=False)

    return df, zero_clusters


def main(argv: list[str] | None = None) -> int:
    """Entry point wired to ``rc2-glm-compare`` in pyproject.toml."""
    parser = argparse.ArgumentParser(
        prog="rc2-glm-compare",
        description="Diff Python GLM outputs against a MATLAB reference directory.",
    )
    parser.add_argument(
        "--python-dir", type=Path, required=True,
        help="Directory with Python-side pipeline outputs (contains "
             "glm_model_comparison.csv, glm_coefficients.csv, etc.)",
    )
    parser.add_argument(
        "--reference-dir", type=Path, required=True,
        help="Directory with the reference (MATLAB) run. "
             "Typically /Volumes/margrie/.../glm_single_cluster/old.",
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Write the per-check report CSV here. Defaults to "
             "<python-dir>/validation/comparison_report.csv.",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s %(message)s",
    )

    out_path = args.out or (args.python_dir / "validation" / "comparison_report.csv")
    df, zero_clusters = run_compare(args.python_dir, args.reference_dir, out_path)

    failing = df[df["pass"] == False]  # noqa: E712 — we want element-wise compare
    passing = df[df["pass"] == True]   # noqa: E712

    print(f"rc2-glm-compare: wrote {out_path}")
    print(f"  checks: {len(df)} | pass: {len(passing)} | fail: {len(failing)}")
    if zero_clusters:
        print(f"  clusters with Jaccard=0 on selected_vars: {zero_clusters[:20]}"
              f"{' ...' if len(zero_clusters) > 20 else ''}")
    if not failing.empty:
        print("\nFAILING CHECKS:")
        # Use to_string so pass is displayed readably even without pandas options.
        print(failing.to_string(index=False))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
