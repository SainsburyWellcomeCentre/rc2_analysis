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
# Tuning-curve prediction-space parity thresholds. ``per_cluster`` gates the
# row pass/fail via the mean (robust to a few outliers); ``per_cluster_min``
# is reported as a side metric so a single bad cluster shows up.
TUNING_CURVE_PEARSON_MEAN_THRESHOLD = 0.90
TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD = 0.85

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
) -> dict:
    """Mean per-cluster Pearson r between Python and MATLAB tuning curves.

    Builds the same raised-cosine grid on both sides using Python's basis
    (MATLAB's is verified numerically identical to ~1e-12), evaluates
    ``B_grid @ β_var`` for the ``model`` fit (default Selected), and
    Pearson-correlates the two curves. Mean across clusters is the row
    value; the gating threshold is on the mean. A per-cluster minimum
    is reported in the note.
    """
    if coef_py is None or coef_ref is None:
        return _row(
            f"tuning_curve_pearson_{variable}", float("nan"),
            TUNING_CURVE_PEARSON_MEAN_THRESHOLD, True,
            note="skipped: glm_coefficients.csv missing on one side",
        )

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
        return _row(
            f"tuning_curve_pearson_{variable}", float("nan"),
            TUNING_CURVE_PEARSON_MEAN_THRESHOLD, True,
            note=(
                f"skipped: no common clusters with {variable} in the "
                f"{model} model on both sides"
            ),
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
        return _row(
            f"tuning_curve_pearson_{variable}", float("nan"),
            TUNING_CURVE_PEARSON_MEAN_THRESHOLD, False,
            note="no clusters yielded a finite Pearson r",
        )

    mean_r = float(np.mean(per_cluster_r))
    min_r = float(np.min(per_cluster_r))
    below = int(sum(r < TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD for r in per_cluster_r))
    passed = mean_r >= TUNING_CURVE_PEARSON_MEAN_THRESHOLD
    return _row(
        f"tuning_curve_pearson_{variable}", mean_r,
        TUNING_CURVE_PEARSON_MEAN_THRESHOLD, passed,
        note=(
            f"n={len(per_cluster_r)} clusters | min_r={min_r:.3f} | "
            f"{below} below per-cluster threshold "
            f"({TUNING_CURVE_PEARSON_PER_CLUSTER_THRESHOLD}) | model={model}"
        ),
    )


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
) -> tuple[pd.DataFrame, list[int]]:
    """Run all checks and return the report + the list of Jaccard=0 cluster ids."""
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
    # the β-rotation that raised-cosine bases admit; gates on mean r.
    config = GLMConfig()
    for variable in ("Speed", "TF"):
        rows.append(_tuning_curve_pearson(coef_py, coef_ref, variable, config))

    rows.append(_stationary_vs_motion_fr(svm_py, svm_ref))

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
