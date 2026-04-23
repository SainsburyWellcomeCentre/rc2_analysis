"""RC2 GLM — Python vs MATLAB parity notebook.

A transparent, audit-ready walk through every step the Python pipeline
takes, with side-by-side checks against the MATLAB reference. The goal
is explicit: be able to trust the Python GLM as the canonical analysis
and stop consulting MATLAB for verification.

Open this file as a notebook in VSCode / Jupyter (it uses the jupytext
``percent`` format: ``# %%`` separates cells, ``# %% [markdown]`` marks
markdown cells). It also runs as a plain script: ``python parity_report.py``
regenerates the comparison report and writes a summary to stdout.

Read the associated design doc: ``summary/rc2_analysis-python-matlab-parity.md``.
Everything in this notebook is traceable to a row in that document.
"""

# %% [markdown]
# # RC2 GLM — Python ↔ MATLAB parity report
#
# **Purpose.** Establish trust in the Python GLM as the canonical analysis
# for the rc2 motion-clouds experiment. Every filter, every threshold,
# every numerical choice is traced to a MATLAB origin and either matched,
# documented as a convention difference, or flagged as a known bug under
# investigation.
#
# **Oracle.** MATLAB reference run on Margrie-lab ceph at
# ``/Volumes/margrie/laura/data transfer for laura/glm_single_cluster/old/``.
# 4 probes: CAA-1123243/244/466/467_rec1.
#
# **Scope of this notebook.** Purely diagnostic: loads the CSVs produced by
# ``rc2-glm`` and the MATLAB reference, runs the parity checks, renders the
# cross-probe figures. It does not re-fit — the aggregated Python run
# lives under ``local_data/motion_clouds/figures/glm_out_all/``.
#
# **What "full trust" means in this document.** The pipeline is trustable
# when every row in the per-probe ``comparison_report.csv`` either (a) passes
# against the MATLAB oracle, or (b) has a named convention difference with a
# rationale you can read from this notebook without re-reading MATLAB source.

# %%
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("parity_report")

# %% [markdown]
# ## 1. Machine + path setup
#
# The paths differ between the two laptops — we detect from ``$HOME``.
# All MATLAB reference CSVs live on Margrie ceph (`/Volumes/margrie/...`);
# the Python aggregated outputs live under ``local_data/motion_clouds/``.

# %%
@dataclass(frozen=True)
class ParityPaths:
    """Resolved paths for this machine.

    ``py_run`` is the aggregated 4-probe directory (concatenated CSVs from
    per-probe runs). ``ml_root`` is the ceph reference directory; missing
    when ceph is not mounted, in which case the notebook can still render
    Python-only sections but parity panels will be skipped.
    """
    home: Path
    py_run: Path
    ml_root: Path
    ml_csv_dir: Path
    formatted_data_dir: Path

    @property
    def ceph_mounted(self) -> bool:
        return self.ml_root.is_dir()


def detect_paths() -> ParityPaths:
    home = Path.home()
    if str(home).startswith("/Users/lauraporta"):
        py_run = Path("/Users/lauraporta/local_data/motion_clouds/figures/glm_out_all")
        formatted = Path("/Users/lauraporta/local_data/motion_clouds/formatted_data")
    elif str(home).startswith("/Users/laura"):
        py_run = Path("/Users/laura/local_data/motion_clouds/figures/glm_out_all")
        formatted = Path("/Users/laura/local_data/motion_clouds/formatted_data")
    else:
        raise RuntimeError(f"unknown machine: home={home}")
    ml_root = Path("/Volumes/margrie/laura/data transfer for laura/glm_single_cluster/old")
    ml_csv = ml_root / "csv"
    return ParityPaths(home=home, py_run=py_run, ml_root=ml_root,
                       ml_csv_dir=ml_csv, formatted_data_dir=formatted)


paths = detect_paths()
print(f"Machine: {paths.home}")
print(f"Python run dir: {paths.py_run}  (exists: {paths.py_run.is_dir()})")
print(f"MATLAB reference: {paths.ml_root}  (mounted: {paths.ceph_mounted})")
print(f"Formatted .mat dir: {paths.formatted_data_dir}  (exists: {paths.formatted_data_dir.is_dir()})")

# %% [markdown]
# ## 2. Experimental paradigm — protocols and velocity channels
#
# The rc2 motion-clouds task runs trials under one of four *protocols*:
#
# | Protocol | Mouse behaviour | Visual flow | Motion lives in |
# |---|---|---|---|
# | **Coupled** | actively locomoting | coupled to running | `filtered_teensy` |
# | **EncoderOnly** | actively locomoting | static / off | `filtered_teensy` |
# | **StageOnly** | passive on moving stage | off | `stage` (stage velocity) |
# | **ReplayOnly** | passive, stage still | replayed from prior trial | `multiplexer_output` ⚠️ |
#
# ⚠️ **ReplayOnly is the sharp-edge case.** The mouse does not translate
# (stage stationary + no running), so the treadmill velocity channel is
# near zero. The ``multiplexer_output`` channel records the *visual-flow*
# velocity that was replayed. Prompt 02.5 fixed the reader so the Python
# pipeline picks the right channel per protocol. Whether the downstream
# motion-*mask* (applied after the velocity trace) also needs this same
# per-protocol awareness is an open question — see prompt 02.8.
#
# We log the protocol distribution per probe from the pipeline's own
# diagnostic CSV.

# %%
def load_trial_channel_report(paths: ParityPaths) -> pd.DataFrame:
    """Concatenate the per-probe trial_channel_report.csv files.

    Writes one row per trial with ``probe_id, trial_id, protocol,
    chosen_channel, n_samples, motion_fraction``. Produced by
    ``rc2_glm.pipeline._emit_trial_channel_report`` during every run.
    """
    runs_root = paths.py_run / "_runs"
    if not runs_root.is_dir():
        raise FileNotFoundError(runs_root)
    dfs = []
    for probe_dir in sorted(runs_root.glob("CAA-*")):
        f = probe_dir / "diagnostics" / "trial_channel_report.csv"
        if f.is_file():
            one = pd.read_csv(f)
            # The per-probe diagnostic doesn't carry probe_id — inject it
            # from the directory name so cross-probe slicing works.
            one.insert(0, "probe_id", probe_dir.name)
            dfs.append(one)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


tcr = load_trial_channel_report(paths)
print(f"Trial-channel rows: {len(tcr)}")
if not tcr.empty:
    print()
    print("Protocol × probe counts:")
    print(tcr.pivot_table(
        index="protocol", columns="probe_id", values="trial_id",
        aggfunc="count", fill_value=0,
    ))
    print()
    print("Chosen velocity channel per protocol (collapsed across probes):")
    print(tcr.groupby(["protocol", "chosen_channel"]).size().unstack(fill_value=0))

# %% [markdown]
# ## 3. Prefilter — which clusters get a GLM
#
# The MATLAB reference runs a paired Wilcoxon test per cluster on
# stationary-vs-motion firing rate for each condition (T, V, VT). The
# 8-branch decision tree categorises clusters as one of:
#
# - ``none_significant``                 — not tuned
# - ``T_only_significant``               — only vestibular
# - ``V_only_significant``               — only visual
# - ``VT_only_significant``              — only multi-sensory
# - ``T_and_V_significant``              — vestibular + visual (additive)
# - ``T_and_VT_significant``             — vestibular + interaction
# - ``V_and_VT_significant``             — visual + interaction
# - ``all_three_significant``            — all three
#
# ``should_run_glm=1`` iff the category is anything other than
# ``none_significant``. Python's prefilter is a line-by-line port of this
# decision tree (``rc2_glm.prefilter``) operating on MATLAB's precomputed
# stationary/motion firing-rate CSV — same inputs, same statistic, same
# categorisation. Therefore we expect **100% category agreement on common
# clusters**.

# %%
def load_prefilter(paths: ParityPaths) -> tuple[pd.DataFrame, pd.DataFrame]:
    py = pd.read_csv(paths.py_run / "prefilter_decision_tree.csv")
    ml = (pd.read_csv(paths.ml_root / "prefilter_decision_tree.csv")
          if paths.ceph_mounted else pd.DataFrame())
    return py, ml


pre_py, pre_ml = load_prefilter(paths)
print(f"Python prefilter rows: {len(pre_py)}  (across {pre_py.probe_id.nunique()} probes)")
if not pre_ml.empty:
    print(f"MATLAB prefilter rows: {len(pre_ml)}  (across {pre_ml.probe_id.nunique()} probes)")
    print()

    # Diagnostic 1: row-count asymmetry. MATLAB only ships clusters that
    # passed its upstream "speed-tuned" filter; Python ships all VISp
    # clusters. This is a scope difference, not a correctness one.
    print("Row counts per probe (Python vs MATLAB prefilter rows):")
    py_counts = pre_py.groupby("probe_id").size().rename("python")
    ml_counts = pre_ml.groupby("probe_id").size().rename("matlab")
    cnts = pd.concat([py_counts, ml_counts], axis=1).fillna(0).astype(int)
    print(cnts)
    print()
    print("  ^ Python = all VISp clusters; MATLAB = clusters that passed its")
    print("    upstream speed-tuning filter. Scope difference, not a bug.")

# %% [markdown]
# ### 3a. Prefilter agreement on common clusters
#
# The confusion matrix below is the MATLAB-vs-Python category assignment
# on the clusters each side retained. Diagonal = agreement, off-diagonal =
# disagreement. The compare tool writes the per-probe + aggregate version
# as ``validation/classification_differences.pdf``; here we compute it
# inline for visibility.

# %%
def prefilter_confusion(pre_py: pd.DataFrame, pre_ml: pd.DataFrame) -> pd.DataFrame:
    merged = pre_py[["probe_id", "cluster_id", "category"]].merge(
        pre_ml[["probe_id", "cluster_id", "category"]],
        on=["probe_id", "cluster_id"], suffixes=("_py", "_ml"),
    )
    print(f"Common (probe_id, cluster_id) pairs: {len(merged)}")
    agree = (merged["category_py"] == merged["category_ml"]).mean()
    print(f"Agreement: {100 * agree:.1f}% ({int(agree * len(merged))}/{len(merged)})")
    disagreements = merged[merged["category_py"] != merged["category_ml"]]
    if len(disagreements):
        print()
        print("DISAGREEMENTS:")
        print(disagreements)
    return merged


if not pre_ml.empty:
    merged_pre = prefilter_confusion(pre_py, pre_ml)

# %% [markdown]
# ## 4. Time-binning and motion mask
#
# Each trial's signals are re-sampled into discrete time bins. Two choices
# here are consequential:
#
# | Choice | Python | MATLAB (reference, 2026-03-11) |
# |---|---|---|
# | bin width | **100 ms** (``config.time_bin_width``) | **20 ms** (diary line 5) |
# | motion mask | ``treadmill_motion_mask`` on the protocol-chosen velocity | ditto, per lab convention |
#
# **Bin width is a known convention difference.** It has no effect on
# forward-selection classification (Spearman invariance) or on tuning-curve
# shape (prediction invariance), but it scales the absolute Poisson
# CV-bits-per-spike by a constant: Python is consistently ~+5.3 bps higher
# than MATLAB. See `summary/rc2_analysis-python-matlab-parity.md`
# (section "CV-bps offset — resolved"). The compare tool reports this as
# an informational row (`null_cv_bps_mean_offset`, `selected_cv_bps_mean_offset`).
#
# **Motion-mask caveat.** On ReplayOnly protocol, there is no treadmill
# translation; motion lives in the visual-flow signal. Whether the mask
# correctly captures ReplayOnly motion bins — and whether MATLAB's
# FullInteraction fit on these clusters works because MATLAB does something
# different — is a pending investigation (prompt 02.8).

# %%
def summarise_bin_counts(paths: ParityPaths) -> pd.DataFrame:
    py_cmp = pd.read_csv(paths.py_run / "glm_model_comparison.csv")
    bin_col = "time_n_bins" if "time_n_bins" in py_cmp.columns else "n_bins"
    return (
        py_cmp.groupby("probe_id")[bin_col]
        .agg(["count", "min", "median", "max"])
        .rename(columns={"count": "n_clusters"})
    )


print("Python per-probe bin-count summary (selected rows):")
print(summarise_bin_counts(paths))

# %% [markdown]
# ## 5. Raised-cosine basis functions
#
# Each continuous predictor (Speed, TF, motion onset time) is expanded into
# raised-cosine basis functions so the GLM can capture non-linear tuning
# without committing to a parametric shape. Speed and TF use the Weber-law
# log-shifted family from MATLAB ``make_raised_cosine_basis``; onset uses
# a causal raised-cosine family.
#
# | Parameter | Speed | TF | Onset |
# |---|---|---|---|
# | n_bases | 5 | 5 | 6 |
# | range | 0 – 50 cm/s | 0 – 7.3 Hz | 0 – 2.0 s |
# | log shift ε | 0.5 | 0.5 | n/a |
# | width | ``delta * 1.5`` | same | same |
#
# Python's implementation (``rc2_glm.basis.raised_cosine_basis``) has a
# docstring assertion that it matches MATLAB output to ~1e-12. This is
# the key invariant — if this ever breaks, tuning-curve Pearson will drop
# across the board and the issue is the basis, not the fit.

# %%
from rc2_glm.basis import onset_kernel_basis, raised_cosine_basis  # noqa: E402
from rc2_glm.config import GLMConfig  # noqa: E402

config = GLMConfig()
speed_grid = np.linspace(config.speed_range[0], config.speed_range[1], 200)
tf_grid = np.linspace(config.tf_range[0], config.tf_range[1], 200)
onset_grid = np.linspace(0.0, config.onset_range[1], 200)

B_speed = raised_cosine_basis(speed_grid, config.n_speed_bases, *config.speed_range)
B_tf = raised_cosine_basis(tf_grid, config.n_tf_bases, *config.tf_range)
B_onset = onset_kernel_basis(onset_grid, config.n_onset_bases, config.onset_range[1])

fig, axes = plt.subplots(1, 3, figsize=(14, 3.2))
for i in range(config.n_speed_bases):
    axes[0].plot(speed_grid, B_speed[:, i], label=f"b{i + 1}")
axes[0].set_xlabel("Speed (cm/s)")
axes[0].set_title(f"Speed bases (n={config.n_speed_bases})")
axes[0].legend(fontsize=7)
for i in range(config.n_tf_bases):
    axes[1].plot(tf_grid, B_tf[:, i], label=f"b{i + 1}")
axes[1].set_xlabel("TF (Hz)")
axes[1].set_title(f"TF bases (n={config.n_tf_bases})")
axes[1].legend(fontsize=7)
for i in range(config.n_onset_bases):
    axes[2].plot(onset_grid, B_onset[:, i], label=f"b{i + 1}")
axes[2].set_xlabel("Time since motion onset (s)")
axes[2].set_title(f"Onset kernel (n={config.n_onset_bases})")
axes[2].legend(fontsize=7)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Forward selection
#
# Hardcastle-style forward-selection builds up the model variable-by-variable:
# at each round, try adding each unselected variable to the current model,
# fit via IRLS on 4/5 of the trials, score on the held-out 5th fold, pick
# the variable with the largest CV-bps improvement, accept if improvement
# ≥ ``delta_bps_threshold`` (0.005 bps). Terminate when no candidate clears
# the threshold.
#
# | Parameter | Value |
# |---|---|
# | n_folds | 5 (trial-level, stratified by condition) |
# | delta_bps_threshold | 0.005 |
# | candidates (main) | Speed, TF, SF, OR |
# | candidates (interactions) | Speed×TF, Speed×SF, Speed×OR, TF×SF, TF×OR, SF×OR |
# | lambda_ridge (main + Selected + Additive) | 0 |
# | lambda_ridge (FullInteraction) | 1.0 |
#
# The Python pipeline uses a ridge of ``1e-3`` for Selected/Additive (commit
# ``11226a8``), which stabilises correlated-basis rotation without biasing
# predictions meaningfully. Without it, Speed/TF β rotate freely on
# correlated raised-cosine bases and raw-β sign parity becomes noise. See
# also ``summary/rc2_analysis-python-matlab-parity.md`` (section
# "Speed/TF parity — resolved").

# %%
def load_model_comparison(paths: ParityPaths) -> tuple[pd.DataFrame, pd.DataFrame]:
    py = pd.read_csv(paths.py_run / "glm_model_comparison.csv")
    py.columns = [c[5:] if c.startswith("time_") else c for c in py.columns]
    ml = pd.DataFrame()
    if paths.ceph_mounted:
        ml = pd.read_csv(paths.ml_csv_dir / "glm_model_comparison.csv")
        ml.columns = [c[5:] if c.startswith("time_") else c for c in ml.columns]
    return py, ml


cmp_py, cmp_ml = load_model_comparison(paths)
print(f"Python model_comparison rows: {len(cmp_py)}  "
      f"(probes: {sorted(cmp_py.probe_id.unique())})")
if not cmp_ml.empty:
    print(f"MATLAB model_comparison rows: {len(cmp_ml)}")
    print()
    # Median / mean number of selected variables per cluster.
    print("Median # selected variables per cluster:")
    print(f"  Python: {cmp_py.n_selected_vars.median():.1f}  "
          f"(max {cmp_py.n_selected_vars.max()})")
    print(f"  MATLAB: {cmp_ml.n_selected_vars.median():.1f}  "
          f"(max {cmp_ml.n_selected_vars.max()})")

# %% [markdown]
# ### 6a. Selected-variables Jaccard per cluster
#
# For each common cluster we compute the Jaccard overlap between the
# Python-selected set and the MATLAB-selected set. Jaccard = 1 is identical
# selection; Jaccard = 0 is entirely disjoint sets; anything in between is
# partial agreement. A histogram of per-cluster Jaccard reveals the shape
# of disagreement.

# %%
def _as_set(s) -> set[str]:
    if not isinstance(s, str):
        return set()
    norm = s.replace(", ", "+").replace(",", "+")
    return {t.strip() for t in norm.split("+") if t.strip() and t.strip() != "Null"}


def jaccard_per_cluster(cmp_py: pd.DataFrame, cmp_ml: pd.DataFrame) -> pd.DataFrame:
    m = cmp_py[["probe_id", "cluster_id", "selected_vars"]].merge(
        cmp_ml[["probe_id", "cluster_id", "selected_vars"]],
        on=["probe_id", "cluster_id"], suffixes=("_py", "_ml"),
    )
    m["py_set"] = m["selected_vars_py"].apply(_as_set)
    m["ml_set"] = m["selected_vars_ml"].apply(_as_set)
    def _j(row):
        a, b = row["py_set"], row["ml_set"]
        if not a and not b:
            return 1.0
        if not (a | b):
            return 1.0
        return len(a & b) / len(a | b)
    m["jaccard"] = m.apply(_j, axis=1)
    return m


if not cmp_ml.empty:
    jac = jaccard_per_cluster(cmp_py, cmp_ml)
    print(f"Jaccard — n={len(jac)} common clusters")
    print(f"  mean   = {jac.jaccard.mean():.3f}")
    print(f"  median = {jac.jaccard.median():.3f}")
    print(f"  identical (jac=1): {int((jac.jaccard == 1.0).sum())}")
    print(f"  partial  (0<jac<1): {int(((jac.jaccard > 0) & (jac.jaccard < 1)).sum())}")
    print(f"  disjoint (jac=0): {int((jac.jaccard == 0).sum())}")

    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(jac.jaccard.values, bins=np.linspace(0, 1, 21),
            edgecolor="black", color="#888")
    ax.set_xlabel("Jaccard(Python, MATLAB) selected_vars")
    ax.set_ylabel("# clusters")
    ax.set_title(f"Per-cluster selection Jaccard (n={len(jac)})")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    plt.show()

    # List clusters with non-trivial disagreement for the record.
    if (jac.jaccard < 1).any():
        print()
        print("Clusters with partial selection disagreement:")
        print(jac[jac.jaccard < 1][[
            "probe_id", "cluster_id", "selected_vars_py", "selected_vars_ml", "jaccard",
        ]].to_string(index=False))

# %% [markdown]
# ## 7. Cross-validated bits-per-spike
#
# For each fitted model we compute the held-out log-likelihood per fold and
# convert to bits-per-spike via
#
# ```
# cv_bps = sum(ll_test) / sum(y_test) / log(2)
# ll_test = sum( y * log(μ) − μ − gammaln(y + 1) )
# ```
#
# Units-wise this is bits of information the model provides about the next
# spike, averaged over test bins. **The absolute value depends on the bin
# width**: at 100 ms (Python) the average per-bin `μ` is larger than at
# 20 ms (MATLAB), so the per-bin Poisson entropy term shifts. This is the
# ~+5.3 bps systematic offset documented under CV-bps in the summary file.
# Spearman rank correlation across clusters is invariant to the offset and
# stays above 0.97.

# %%
if not cmp_ml.empty:
    m = cmp_py[["probe_id", "cluster_id", "Null_cv_bps", "Selected_cv_bps"]].merge(
        cmp_ml[["probe_id", "cluster_id", "Null_cv_bps", "Selected_cv_bps"]],
        on=["probe_id", "cluster_id"], suffixes=("_py", "_ml"),
    ).dropna()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, col in zip(axes, ("Null", "Selected")):
        x = m[f"{col}_cv_bps_ml"].values
        y = m[f"{col}_cv_bps_py"].values
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], "k:", linewidth=1, alpha=0.4, label="y=x")
        offset = float((y - x).mean())
        ax.plot([lo, hi], [lo + offset, hi + offset], color="#c33", linewidth=0.8,
                label=f"y=x+{offset:.2f}")
        ax.scatter(x, y, s=12, alpha=0.7, color="#334")
        ax.set_xlabel(f"MATLAB {col} cv_bps")
        ax.set_ylabel(f"Python {col} cv_bps")
        ax.set_title(f"{col}: offset={offset:+.3f} bps (bin-width convention)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()

    # Spearman is invariant to the offset — this is the parity gate.
    from scipy.stats import spearmanr
    for col in ("Null", "Selected"):
        x = m[f"{col}_cv_bps_ml"].values
        y = m[f"{col}_cv_bps_py"].values
        rho, _ = spearmanr(x, y)
        print(f"  Spearman {col}_cv_bps: {rho:.3f}  (n={len(x)})")

# %% [markdown]
# ## 8. Tuning-curve parity (prediction-space)
#
# Raw-β sign agreement on the Speed/TF raised-cosine bases is not
# identifiable — two β vectors differing by a rotation can produce the
# same predicted tuning curve. Prompt 02.7 moved the gate from
# ``sign(β_speed)`` to Pearson correlation of the predicted tuning curves
# on a shared grid. Gate was further refined in commit ``e677b9f``:
# ``fraction(per_cluster_r ≥ 0.85) ≥ 0.80`` so single-cluster outliers
# at small n (n=6 TF-selected) don't trip a row spuriously.
#
# **Known caveat (future work).** The tuning curves are evaluated at
# ``t_since_onset = 1.5 s`` (steady-state of the onset kernel) and with
# the other continuous variable (TF when plotting Speed, and vice versa)
# held at its empirical mean. Both are conventions, not facts about the
# cells. Prompt 02.8 (future) adds an opt-in "trial-averaged" mode that
# marginalises the onset basis and the other continuous variables over
# their empirical per-cluster distributions. Until then, treat tuning
# curves as **steady-state**, not as trial-averaged.

# %%
from rc2_glm.compare import _tuning_curve_pearson, _variable_betas_per_cluster  # noqa: E402


def tuning_pearson_per_cluster(
    coef_py: pd.DataFrame, coef_ml: pd.DataFrame, variable: str,
) -> pd.DataFrame:
    """Per-cluster Pearson r — exposed so the notebook can plot the
    distribution rather than the single row value in comparison_report.csv."""
    rows_gate, rows_info = _tuning_curve_pearson(coef_py, coef_ml, variable, config)
    # _tuning_curve_pearson hides the per-cluster list; re-derive here so
    # we can histogram and identify outliers.
    n_bases = config.n_speed_bases if variable == "Speed" else config.n_tf_bases
    lo, hi = config.speed_range if variable == "Speed" else config.tf_range
    py_map = _variable_betas_per_cluster(coef_py, variable, "Selected", n_bases)
    ml_map = _variable_betas_per_cluster(coef_ml, variable, "Selected", n_bases)
    x = np.linspace(float(lo), float(hi), 100)
    B = raised_cosine_basis(x, n_bases, float(lo), float(hi))
    rows = []
    for key in sorted(set(py_map) & set(ml_map)):
        c_py = B @ py_map[key]
        c_ml = B @ ml_map[key]
        if np.std(c_py) < 1e-12 or np.std(c_ml) < 1e-12:
            continue
        r = np.corrcoef(c_py, c_ml)[0, 1]
        if np.isfinite(r):
            rows.append({
                "probe_id": key[0], "cluster_id": key[1], "variable": variable,
                "pearson_r": float(r),
            })
    return pd.DataFrame(rows)


def load_coefficients(paths: ParityPaths) -> tuple[pd.DataFrame, pd.DataFrame]:
    py = pd.read_csv(paths.py_run / "glm_coefficients.csv")
    ml = (pd.read_csv(paths.ml_csv_dir / "glm_coefficients.csv")
          if paths.ceph_mounted else pd.DataFrame())
    return py, ml


coef_py, coef_ml = load_coefficients(paths)
print(f"Python coefficients rows: {len(coef_py)}")
if not coef_ml.empty:
    print(f"MATLAB coefficients rows: {len(coef_ml)}")
    print()
    for var in ("Speed", "TF"):
        per_cluster = tuning_pearson_per_cluster(coef_py, coef_ml, var)
        n = len(per_cluster)
        r = per_cluster.pearson_r.values
        frac = float((r >= 0.85).mean()) if n else float("nan")
        print(f"  {var}: n={n} | mean r={r.mean():.3f} | median={np.median(r):.3f} | "
              f"frac(r>=0.85)={frac:.3f} | min={r.min():.3f}")

# %%
if not coef_ml.empty:
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
    for ax, var in zip(axes, ("Speed", "TF")):
        per_cluster = tuning_pearson_per_cluster(coef_py, coef_ml, var)
        if per_cluster.empty:
            ax.text(0.5, 0.5, f"no common {var}", transform=ax.transAxes, ha="center")
            continue
        r = per_cluster.pearson_r.values
        ax.hist(r, bins=np.linspace(-1, 1, 41), edgecolor="black", color="#668")
        ax.axvline(0.85, color="#c33", linestyle="--", label="per-cluster bar 0.85")
        ax.set_xlabel("Pearson r (Python, MATLAB tuning curve)")
        ax.set_ylabel("# clusters")
        ax.set_title(f"{var}: n={len(r)}, frac≥0.85 = {(r >= 0.85).mean():.3f}")
        ax.legend(fontsize=8)
    fig.tight_layout()
    plt.show()

# %% [markdown]
# ## 9. Cross-probe classification figures (recap)
#
# The compare tool writes two headline figures under
# ``{py_run}/validation/``. Re-render them inline here so this notebook
# is the single document a reviewer needs.

# %%
from rc2_glm.compare import (  # noqa: E402
    _plot_classification_differences,
    _plot_selected_vars_differences,
)

validation_dir = paths.py_run / "validation"
validation_dir.mkdir(exist_ok=True)

if not pre_ml.empty:
    _plot_classification_differences(
        pre_py, pre_ml,
        out_pdf=validation_dir / "classification_differences.pdf",
        out_disagreements_csv=validation_dir / "classification_disagreements.csv",
    )
    print(f"wrote {validation_dir / 'classification_differences.pdf'}")

if not cmp_ml.empty:
    cmp_py_raw = pd.read_csv(paths.py_run / "glm_model_comparison.csv")
    cmp_ml_raw = pd.read_csv(paths.ml_csv_dir / "glm_model_comparison.csv")
    _plot_selected_vars_differences(
        cmp_py_raw, cmp_ml_raw,
        out_pdf=validation_dir / "selected_vars_differences.pdf",
        out_disagreements_csv=validation_dir / "selected_vars_disagreements.csv",
    )
    print(f"wrote {validation_dir / 'selected_vars_differences.pdf'}")

# %% [markdown]
# ## 10. Known divergences (audit table)
#
# | Topic | Divergence | Resolution | See |
# |---|---|---|---|
# | Bin width | Python 100 ms vs MATLAB 20 ms | Document, don't force — CV-bps mean offset is informational | summary §CV-bps |
# | SF / OR coefficient naming | MATLAB uses integer indices for FullInteraction | Python uses level values | summary §"Coefficient naming" |
# | Null + FullInteraction β export | Python CSV only exports Selected + Additive | Minor; Selected covers the scientific result | summary §"Exported models" |
# | Upstream speed-tuning filter | MATLAB prefilter restricts to speed-tuned clusters upstream | Python keeps all VISp (broader scope) | this notebook §3 |
# | Forward-selection boundary | Python picks ~5 more interaction terms near Δ-bps=0.005 | Ridge + RNG; near-boundary flips | summary §"Jaccard" |
# | FullInteraction skip rate | Python silently drops clusters with p ≥ n | Hypothesis: ReplayOnly motion mask | prompt 02.8 |
# | Tuning-curve steady-state | Currently evaluated at t=1.5s post-onset | Trial-averaged mode planned | prompt 02.8/8.6 |
#
# **Everything above is either benign (convention), compensated (ridge
# stabilisation, prediction-space metric), or queued as a named
# investigation with a reviewer prompt.**

# %% [markdown]
# ## 11. Trust checklist — can we move on?
#
# For Laura's "fully trust Python" milestone:
#
# 4-probe aggregate (101 common prefilter rows, 88 common GLM rows,
# 18/18 parity checks pass):
#
# - [x] Prefilter category agreement ≥ 90% — **100%** (101/101)
# - [x] Forward-selection Jaccard ≥ 0.70 — **0.921** mean, 0 disjoint, 13 partial
# - [x] CV-bps rank correlation ≥ 0.85 — Null **0.977**, Selected **0.979**
# - [x] Selected−Null direction agreement ≥ 0.95 — **1.000**
# - [x] Categorical coefficient-sign ≥ 0.85 — SF **0.938**, OR **1.000**
# - [x] Tuning-curve frac(r ≥ 0.85) ≥ 0.80 — Speed **0.803**, TF **0.850**
# - [x] Cross-probe parity figures (4-probe confusion + variable bars) present
# - [ ] ReplayOnly motion-mask diagnostic (prompt 02.8 — queued, not gating)
# - [ ] Trial-averaged tuning-curve mode (prompt 02.8/8.6 — queued, cosmetic)
#
# Two open items remain, both documented as follow-up prompts. The
# analysis pipeline is otherwise trustable — the divergences from MATLAB
# are either convention differences (bin width), identifiability
# artefacts already handled in the parity metric (raw-β sign →
# tuning-curve Pearson), or per-cluster near-boundary flips that the
# parity summary lists by cluster ID.

# %% [markdown]
# ## 12. Regenerate this report
#
# The `rc2-glm-compare` CLI regenerates the CSV + figures in one shot:
#
# ```bash
# rc2-glm-compare \
#   --python-dir  /Users/lauraporta/local_data/motion_clouds/figures/glm_out_all \
#   --reference-dir "/Volumes/margrie/laura/data transfer for laura/glm_single_cluster/old" \
#   --out /Users/lauraporta/local_data/motion_clouds/figures/glm_out_all/validation/comparison_report.csv
# ```
#
# The `--python-dir` can point at any Python run — a single-probe dir, or
# the aggregated multi-probe dir built by concatenating the per-probe
# CSVs. Comparison rows, the classification confusion matrix, the
# selected-vars bar chart, and the per-cluster disagreements CSV all land
# in `{python-dir}/validation/`.
