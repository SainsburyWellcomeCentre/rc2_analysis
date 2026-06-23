#!/usr/bin/env python3
"""FENS combined per-cluster figure — 3 rows x 5 variable-columns.

For each requested cluster, ONE figure (one PDF/PNG per cluster):

    columns = Speed | Acceleration | TF | OR | SF      (the 5 tuning variables)

    row 1  OBSERVED tuning
        per-trial-per-bin mean +/- SEM over 20 equal-count (5%-quantile) bins,
        one trace per condition, plus the R^2-selected tuning fit (von Mises for
        OR; linear / quadratic / cubic / Gaussian / asym-Gaussian / sigmoid for
        the rest, family by highest R^2-on-the-mean-curve), with R^2(mean) and
        the significance p READ from the cached tuning_fits.csv — this figure
        never bootstraps. TF/SF/OR are in that cache; Speed/Acc are not, so their
        fit is derived by least-squares only (R^2 shown, p as n/a).

    row 2  GLM KERNELS for the variables forward selection chose
        raised-cosine / circular reconstruction (basis @ fitted beta) from the
        run's glm_coefficients.csv, exactly as cluster_<id>_kernels.pdf draws
        them. A column whose main effect was NOT forward-selected is marked
        "not selected". Selected interaction terms are noted in the row label.

    row 3  MODEL-PREDICTED tuning
        reconstruct each trial's firing rate with the Selected model at the
        native 20 ms bin (mu = exp(X @ beta) with each bin's ACTUAL covariates;
        no marginal sweep, nothing held fixed), then bin those reconstructed
        per-bin rates into the SAME 5%-quantile bins as row 1 and draw the
        per-trial-per-bin mean +/- SEM. ONE reconstruction principle across all
        five columns, so row 1 (observed) and row 3 (predicted) are a true
        apples-to-apples comparison; row 2 carries the isolated kernel.

Conditions overlaid per variable (Laura: keep VF / T / VT):
    Speed         -> T_Vstatic + VT   (Speed == 0 in V)
    Acceleration  -> T_Vstatic + VT   (the MATLAB accel cache has no V)
    TF / SF / OR  -> V + VT            (TF == 0, SF/OR == NaN in T_Vstatic)

Faithfulness. Row 3 uses the EXACT production betas read from
glm_coefficients.csv (model == "Selected"), aligned by column name to a design
matrix rebuilt with the production config (make_config_histbase_all). The align
step fails loud if any saved coefficient has no column to land on (a zero-filled
trained term — the standing `grep ZERO-FILLED` correctness check), so the
per-bin reconstruction is byte-faithful to what the cohort fit produced.

Standalone: reuses make_fens_poster_figures + rc2_glm; modifies neither. The run
directory is a single constant (RUN_DIR) — flip it to the next run when it lands.

Usage:
    python make_fens_combined_cluster_figure.py \
        --probe CAA-1124371_rec1_rec2_rec3 --clusters 97 14 62 311
"""
from __future__ import annotations

# Pin BLAS/OMP to a single thread BEFORE numpy imports, matching
# make_fens_poster_figures (keeps any numeric path reproducible).
import os as _os

for _v in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    _os.environ.setdefault(_v, "1")

import argparse
import logging
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rc2_glm.io import load_probe_data
from rc2_glm.precomputed_bins import load_precomputed_bin_edges
from rc2_glm.plots import COND_COLORS, _kernel_for_var, _pooled_quantile_edges
from rc2_glm.rf_sf_or import load_rf_sf_or
from rc2_glm.time_binning import bin_cluster
from rc2_glm.tuning_significance import (
    evaluate as eval_tuning_fit,
    per_trial_bin_matrix,
    rsq_against_mean,
    tuning_significance,
)

# Reuse the FENS helpers verbatim (observed binning, design prep, selection
# read, family set, axis labels) so this figure shares ONE code path with the
# poster figures and the production fits.
from make_fens_poster_figures import (  # noqa: E402
    TUN_LINEAR_FAMILIES,
    TUN_XLABEL,
    _design,
    _load_fits_lookup,
    _observed_value_tuning,
    _sig_from_cache,
    _tuning_stats,
    prepare_cluster_design,
)
from run_glm_goggles_rf_sfor_20ms import (  # noqa: E402
    FORMATTED_DIR,
    _lookup,
    make_config_histbase_all,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("fens_combined")

HOME = Path.home()
_GLM_DIR = HOME / "local_data" / "motion_clouds" / "figures" / "glm"

# --- the ONE run-directory constant -------------------------------------- #
# The histbase "_all" goggles run (forward-selection history + coefficients).
# Flip this line to the next run when it lands.
RUN_DIR = _GLM_DIR / "current_rf_sfor_20ms_histbase_goggles_all"

OUT_DIR = _GLM_DIR / "FENS_figures_poster" / "combined_cluster"

# Cached tuning fits (selected family + params + the bootstrap p) from a prior
# FENS tuning-browse run made with the SAME settings — mean+SEM display, R²-mean
# family selection. The p-value is READ from here; this figure never bootstraps.
# The cache covers TF/SF/OR (V/VT); Speed/Acc are not in it, so their fit is
# derived by least-squares only (R² shown, no p — n_reps=0, no shuffle).
FITS_CSV = (_GLM_DIR / "FENS_figures_poster" / "tuning_browse_mean_sem_selR2mean"
            / "tuning_fits.csv")

# Columns, in the order Laura asked: Speed, Acc, TF, OR, SF.
#   key       -> tuning_significance value + pc cache + edge logic
#   value_col -> the per-bin covariate column in the binned df
#   kvar      -> the variable name _kernel_for_var expects
#   conds     -> conditions to overlay (where the variable is defined)
COLUMNS = [
    dict(key="speed", value_col="speed", kvar="Speed",
         label="Speed", conds=("T_Vstatic", "VT")),
    dict(key="accel", value_col="acceleration", kvar="Acceleration",
         label="Acceleration", conds=("T_Vstatic", "VT")),
    dict(key="tf", value_col="tf", kvar="TF",
         label="TF", conds=("V", "VT")),
    dict(key="or", value_col="orientation", kvar="OR",
         label="Orientation", conds=("V", "VT")),
    dict(key="sf", value_col="sf", kvar="SF",
         label="SF", conds=("V", "VT")),
]


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_probe(probe_id: str, config):
    """Load probe + RF lookup + precomputed (MATLAB) tuning cache once."""
    probe = load_probe_data(
        FORMATTED_DIR / f"{probe_id}.mat", config=config,
        stimulus_lookup=_lookup(), cluster_set="selected",
    )
    rf_lookup = load_rf_sf_or(
        config.rf_sf_or_parquet_dir, probe.probe_id,
        min_concentration=getattr(config, "rf_min_concentration", 0.0),
    )
    by_id = {c.cluster_id: c for c in probe.clusters}
    pc = load_precomputed_bin_edges(FORMATTED_DIR / f"{probe_id}.mat")
    return probe, by_id, rf_lookup, pc


def load_selected_coefficients(probe_id: str, cluster_id: int) -> dict[str, float]:
    """{coefficient_name: estimate} for the Selected model of one cluster,
    read from the run's per-probe glm_coefficients.csv."""
    csv = RUN_DIR / "_runs" / probe_id / "glm_coefficients.csv"
    if not csv.exists():
        raise FileNotFoundError(csv)
    df = pd.read_csv(csv)
    sub = df[(df["cluster_id"] == cluster_id) & (df["model"] == "Selected")]
    if sub.empty:
        raise ValueError(f"no Selected coefficients for {probe_id} cl {cluster_id}")
    return {str(r.coefficient): float(r.estimate) for r in sub.itertuples()}


def load_selection_vars(probe_id: str, cluster_id: int) -> list[str]:
    """Forward-selected terms (selection order) from glm_selection_history.csv
    — the var list `_design` needs to rebuild the Selected design matrix."""
    csv = RUN_DIR / "_runs" / probe_id / "glm_selection_history.csv"
    sh = pd.read_csv(csv)
    sub = sh[sh["cluster_id"] == cluster_id].copy()
    added = sub["added"].astype(str).str.lower() == "true"
    sub = sub[added].sort_values("round")
    return [str(c) for c in sub["best_candidate"].tolist() if str(c) not in ("", "nan")]


# --------------------------------------------------------------------------- #
# Row 3 — reconstruct each trial's per-bin FR with the Selected model
# --------------------------------------------------------------------------- #
def selected_per_bin_prediction(df, config, probe_id, cluster_id):
    """Per-20ms-bin predicted FR (Hz) for the Selected model, using each bin's
    actual covariates and the EXACT saved betas (no refit).

    Returns (pred_fr_all_rows, parity) where parity is a dict describing the
    column-set check (the zero-filled-term guard).
    """
    prep = prepare_cluster_design(df, config)
    selected_vars = load_selection_vars(probe_id, cluster_id)
    X, col_names = _design(prep, selected_vars, config)
    coef = load_selected_coefficients(probe_id, cluster_id)

    # Align saved betas to the rebuilt design BY NAME. A saved coefficient with
    # no column to land on, or a built column with no saved beta, is a parity
    # break — the same failure class the production `_align_prediction_columns`
    # warns about (a silently zero-filled trained term).
    missing_in_design = [n for n in coef if n not in col_names]   # saved -> no column
    missing_in_coef = [n for n in col_names if n not in coef]     # column -> zero-filled
    beta = np.array([coef.get(n, 0.0) for n in col_names], dtype=np.float64)

    # exp(X @ beta) == count/bin / bin_width == rate in Hz (offset = log(bw)).
    pred_fr = np.exp(np.clip(X @ beta, -20.0, 20.0))
    parity = dict(
        n_saved=len(coef), n_design=len(col_names),
        missing_in_design=missing_in_design, zero_filled=missing_in_coef,
        ok=(not missing_in_design and not missing_in_coef),
    )
    return pred_fr, parity, selected_vars, coef


def _accel_edges_from_centres(centres: np.ndarray) -> np.ndarray | None:
    """Reconstruct bin edges from cache bin centres (midpoints, ends mirrored).
    The accel cache exposes centres but not edges; row 1 (observed) uses the
    cache centres, so row 3 must digitise into the matching edges."""
    if centres is None or len(centres) < 2:
        return None
    c = np.asarray(centres, float)
    inner = 0.5 * (c[:-1] + c[1:])
    first = c[0] - (inner[0] - c[0])
    last = c[-1] + (c[-1] - inner[-1])
    return np.concatenate([[first], inner, [last]])


def edges_for(col, cond, df, pc):
    """The SAME 5%-quantile bin edges row 1 uses for (variable, condition):
    MATLAB cache edges for Speed/TF, cache-centre-derived edges for Accel,
    pooled (V+VT) quantile edges for SF/OR."""
    key, value_col = col["key"], col["value_col"]
    if key == "speed":
        return pc.speed_edges(cond)
    if key == "tf":
        return pc.tf_edges(cond)
    if key == "accel":
        return _accel_edges_from_centres(pc.accel_centres(cond))
    # SF / OR: pooled across the visual conditions (V+VT), the project
    # convention so bin k is the same interval in both — matches
    # _observed_value_tuning for these continuous rf_local covariates.
    pooled = df.loc[df["condition"].isin(("V", "VT")), value_col].to_numpy(float)
    pooled = pooled[np.isfinite(pooled)]
    if pooled.size < 20:
        return None
    edges, _ = _pooled_quantile_edges(pooled, n_bins=20)
    return edges


def predicted_tuning_stats(df_pred, value_col, cond, edges, bw):
    """Per-trial-per-bin mean/SEM of the model-predicted FR, binned into `edges`
    — identical binning to the observed row (per_trial_bin_matrix), just the
    predicted rate (carried as spike_count = pred_fr * bw) instead of spikes."""
    if edges is None:
        return None
    matrix, centres = per_trial_bin_matrix(
        df_pred, value_col, bw, condition=cond, n_bins=len(edges) - 1, edges=edges)
    if matrix is None:
        return None
    return _tuning_stats(matrix, centres, source="model_predicted")


# --------------------------------------------------------------------------- #
# Panel renderers
# --------------------------------------------------------------------------- #
def render_row1(ax, df, col, pc, bw, cluster_id, probe_id, fits_lookup):
    """Observed mean+/-SEM per condition + the R^2-selected fit.

    The fit + p are READ from the cached tuning_fits.csv when present (TF/SF/OR);
    otherwise (Speed/Acc) the R^2-selected family is derived by least-squares with
    n_reps=0 — fit only, NO shuffle bootstrap, p shown as n/a.
    """
    key, value_col, conds = col["key"], col["value_col"], col["conds"]
    kind = "circular" if key == "or" else "linear"
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(TUN_XLABEL[key], fontsize=8)
    notes = []
    drew = False
    for cond in conds:
        t = _observed_value_tuning(df, value_col, key, cond, bw,
                                   cluster_id=cluster_id, pc=pc)
        if t is None:
            continue
        cen = t["centres"]
        good = np.isfinite(t["mean"])
        if not good.any():
            continue
        drew = True
        ax.errorbar(cen[good], t["mean"][good], yerr=t["sem"][good],
                    fmt="o-", color=COND_COLORS[cond], ms=3.5, lw=1.1,
                    capsize=2.5, elinewidth=0.8, label=cond, zorder=2)
        # R^2-selected tuning fit — cached (read the saved p, no bootstrap) or,
        # for the uncached Speed/Acc cells, least-squares only (n_reps=0).
        cached = (None if fits_lookup is None
                  else fits_lookup.get((probe_id, cluster_id, key, cond)))
        try:
            if cached is not None:
                sig = _sig_from_cache(cached)
                if not np.isfinite(sig.rsq_mean) and sig.params is not None:
                    sig.rsq_mean = rsq_against_mean(
                        t["matrix"], cen, sig.best_model, sig.params)
            else:
                sig = tuning_significance(
                    t["matrix"], cen, value=key, condition=cond, kind=kind,
                    aggregate="flat", select_criterion="rsq_mean",
                    n_reps=0, linear_families=TUN_LINEAR_FAMILIES)
        except Exception as exc:  # degenerate cell -> skip the curve, keep points
            log.warning("fit failed %s %s cl%s: %s", key, cond, cluster_id, exc)
            sig = None
        if sig is not None and sig.best_model is not None:
            xs = np.linspace(float(np.nanmin(cen)), float(np.nanmax(cen)), 200)
            ax.plot(xs, eval_tuning_fit(
                dict(name=sig.best_model, params=sig.params), xs),
                color=COND_COLORS[cond], lw=1.7, alpha=0.9, zorder=3)
            if np.isfinite(sig.p):
                pstr = "p<0.001" if sig.p < 1e-3 else f"p={sig.p:.3f}"
            else:
                pstr = "p:n/a"
            notes.append((cond, f"{sig.best_model} R²m={sig.rsq_mean:.2f} {pstr}"))
    if not drew:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="#888", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        return
    if notes:
        txt = "\n".join(f"{c}: {s}" for c, s in notes)
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, ha="left", va="top",
                fontsize=5.4, zorder=5,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=0.85))


def render_row2(ax, col, coef, config, interactions, accel_unscale=None):
    """Forward-selected kernel for this column's main effect (or 'not selected').
    `interactions` (mutated) collects selected interaction groups for the label.
    `accel_unscale=(mean, std)` maps the Acceleration kernel's z-score grid back
    to cm/s² so its x-axis matches the observed/predicted Acc rows."""
    kvar = col["kvar"]
    coef_rows = pd.DataFrame(
        {"coefficient": list(coef.keys()), "estimate": list(coef.values())})
    ax.spines[["top", "right"]].set_visible(False)
    out = _kernel_for_var(coef_rows, kvar, config)
    if out is None:
        ax.text(0.5, 0.5, "not selected", ha="center", va="center",
                transform=ax.transAxes, color="#aaa", fontsize=9, style="italic")
        ax.set_xticks([]); ax.set_yticks([])
    else:
        x, y, xlabel = out
        if col["key"] == "accel" and accel_unscale is not None:
            # invert prepare_cluster_design's z-scoring: cm/s² = z*std + mean
            a_mean, a_std = accel_unscale
            x = np.asarray(x, float) * a_std + a_mean
            xlabel = TUN_XLABEL["accel"]
        ax.axhline(0.0, color="0.8", lw=0.7, zorder=1)
        ax.plot(x, y, color="black", lw=1.8, zorder=2)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("kernel (log-rate)", fontsize=7)
    # record any selected interactions that touch this main effect
    for name in coef:
        if "_x_" in name:
            grp = _interaction_group(name)
            if grp:
                interactions.add(grp)


def render_row3(ax, df_pred, col, pc, df, bw):
    """Model-predicted mean+/-SEM tuning per condition (binned like row 1)."""
    key, value_col, conds = col["key"], col["value_col"], col["conds"]
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(TUN_XLABEL[key], fontsize=8)
    drew = False
    for cond in conds:
        edges = edges_for(col, cond, df, pc)
        stats = predicted_tuning_stats(df_pred, value_col, cond, edges, bw)
        if stats is None:
            continue
        cen = stats["centres"]
        good = np.isfinite(stats["mean"])
        if not good.any():
            continue
        drew = True
        ax.errorbar(cen[good], stats["mean"][good], yerr=stats["sem"][good],
                    fmt="o-", color=COND_COLORS[cond], ms=3.5, lw=1.1,
                    capsize=2.5, elinewidth=0.8, label=cond, zorder=2)
    if not drew:
        ax.text(0.5, 0.5, "no prediction", ha="center", va="center",
                transform=ax.transAxes, color="#888", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])


def _interaction_group(name: str) -> str | None:
    """'Spd2_x_SF_1' -> 'Spd x SF' (strip basis indices), for the row-2 label."""
    parts = name.split("_x_")
    if len(parts) != 2:
        return None
    def base(p):
        return re.sub(r"_?\d+$", "", p)
    return f"{base(parts[0])} x {base(parts[1])}"


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def make_cluster_figure(probe, by_id, rf_lookup, pc, config, probe_id,
                        cluster_id, fits_lookup):
    if cluster_id not in by_id:
        log.warning("cluster %s not in selected set for %s — skipping",
                    cluster_id, probe_id)
        return None
    bw = float(config.time_bin_width)
    df = bin_cluster(probe, by_id[cluster_id], rf_lookup=rf_lookup)

    # Accel z-score -> cm/s² inverse (mean/std over finite motion rows, the same
    # statistics prepare_cluster_design uses) for the Acceleration kernel axis.
    acc_raw = df["acceleration"].to_numpy(float)
    afin = np.isfinite(acc_raw) & (df["condition"] != "stationary").to_numpy()
    accel_unscale = ((float(acc_raw[afin].mean()),
                      float(acc_raw[afin].std(ddof=0)) or 1.0)
                     if afin.sum() >= 1 else None)

    # Row-3 reconstruction (saved betas, parity-checked).
    pred_fr, parity, selected_vars, coef = selected_per_bin_prediction(
        df, config, probe_id, cluster_id)
    if not parity["ok"]:
        log.warning("PARITY cl%s: zero_filled=%s missing_in_design=%s",
                    cluster_id, parity["zero_filled"], parity["missing_in_design"])
    else:
        log.info("parity OK cl%s (%d coeffs == %d columns)",
                 cluster_id, parity["n_saved"], parity["n_design"])
    df_pred = df.copy()
    df_pred["spike_count"] = pred_fr * bw  # so per_trial_bin_matrix FR == pred_fr

    fig, axes = plt.subplots(3, len(COLUMNS), figsize=(4.0 * len(COLUMNS), 9.6),
                             constrained_layout=True)
    interactions: set[str] = set()
    for c, col in enumerate(COLUMNS):
        if c == 0:
            axes[0, c].set_ylabel("OBSERVED\nFR (Hz)\nmean +/- SEM", fontsize=8)
            axes[1, c].set_ylabel("GLM KERNEL\n(selected)", fontsize=8)
            axes[2, c].set_ylabel("PREDICTED\nFR (Hz)\nmean +/- SEM", fontsize=8)
        axes[0, c].set_title(col["label"], fontsize=11, fontweight="bold")
        render_row1(axes[0, c], df, col, pc, bw, cluster_id, probe_id, fits_lookup)
        render_row2(axes[1, c], col, coef, config, interactions,
                    accel_unscale=accel_unscale)
        render_row3(axes[2, c], df_pred, col, pc, df, bw)

    # Share y within each tuning column between observed (row0) and predicted
    # (row2) so the comparison is on one scale; kernels (row1) keep their own.
    for c in range(len(COLUMNS)):
        his = [axes[r, c].get_ylim()[1] for r in (0, 2)
               if np.isfinite(axes[r, c].get_ylim()[1])]
        lo = [axes[r, c].get_ylim()[0] for r in (0, 2)
              if np.isfinite(axes[r, c].get_ylim()[0])]
        if his:
            hi = max(his)
            for r in (0, 2):
                axes[r, c].set_ylim(min(0.0, min(lo)) if lo else 0.0, hi * 1.05)

    handles = [plt.Line2D([0], [0], marker="o", color=COND_COLORS[k], lw=1.1,
                          label={"T_Vstatic": "T (vestibular)", "V": "VF (visual)",
                                 "VT": "VF+T"}[k])
               for k in ("T_Vstatic", "V", "VT")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.012))

    sel = "+".join(selected_vars) if selected_vars else "Null"
    int_txt = ("   interactions: " + ", ".join(sorted(interactions))) if interactions else ""
    parity_txt = "" if parity["ok"] else "   [!] PARITY: zero-filled " + ",".join(parity["zero_filled"])
    probe_short = probe_id.split("_rec")[0]
    fig.suptitle(
        f"{probe_short}  cluster {cluster_id}   "
        f"observed tuning  |  forward-selected GLM kernels  |  model-predicted tuning"
        f"\nselected: {sel}{int_txt}   [{RUN_DIR.name}]{parity_txt}",
        fontsize=11, fontweight="bold")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUT_DIR / f"combined_{probe_short}_cl{cluster_id}"
    for ext in ("pdf", "png"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", stem)
    return stem


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--probe", default="CAA-1124371_rec1_rec2_rec3")
    ap.add_argument("--clusters", type=int, nargs="+", default=[97, 14, 62, 311])
    ap.add_argument("--fits-csv", default=str(FITS_CSV),
                    help="Cached tuning_fits.csv to READ fit + p from (no bootstrap).")
    args = ap.parse_args()

    config = make_config_histbase_all()
    fits_lookup = None
    fits_csv = Path(args.fits_csv)
    if fits_csv.exists():
        fits_lookup = _load_fits_lookup(fits_csv)
        log.info("read %d cached tuning fits from %s (no bootstrap)",
                 len(fits_lookup), fits_csv)
    else:
        log.warning("no cached fits at %s — Speed/Acc/TF/OR/SF fits derived by "
                    "least-squares (R² only, no p)", fits_csv)

    log.info("loading %s ...", args.probe)
    probe, by_id, rf_lookup, pc = load_probe(args.probe, config)

    written = []
    for cid in args.clusters:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stem = make_cluster_figure(probe, by_id, rf_lookup, pc, config,
                                       args.probe, cid, fits_lookup)
        if stem is not None:
            written.append(stem)
    log.info("done — %d figure(s) in %s", len(written), OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
