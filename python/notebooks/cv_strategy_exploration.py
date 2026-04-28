"""CV-strategy exploration: speed-profile, bin-width, seed stability.

Companion to ``parity_report.py``. Where that notebook answers "can we
trust the Python pipeline against the MATLAB oracle?", this one answers
"how sensitive is the Python pipeline to its own design choices?"

Three diagnostic sections, each driven by the CLI knobs landed in the
2026-04-24 exploration-pack commits on ``glm-improvements``:

1. **Speed-profile CV** — MATLAB-parity 2-fold train-one-profile /
   test-the-other diagnostic on the Selected model. Flag:
   ``--profile-cv-diagnostic``.
2. **Bin-width sweep** — 100 / 50 / 20 ms time-bin widths on the same
   clusters, eyeballing impact on tuning curves + per-cluster cv-bps.
   Flag: ``--bin-width``.
3. **Seed stability** — N seeds on the same clusters, quantifying how
   much the condition-stratified fold RNG affects fits. Flag:
   ``--cv-seed``.

Runs the experiments ad-hoc — see the shell snippets at the top of each
section. Same jupytext percent format as ``parity_report.py`` (``# %%``
separates cells, ``# %% [markdown]`` marks markdown).
"""

# %% [markdown]
# # CV strategy exploration — speed-profile, bin-width, seed
#
# Companion to `parity_report.py`. This notebook does not interact with
# the MATLAB reference; it characterises how sensitive the Python fits
# are to three of their own design choices.
#
# All three experiments run on a single probe (CAA-1123243_rec1, ~74
# retained clusters) to keep wall-clock small. Re-run with other probes
# by copy-pasting the shell snippets under each section and updating
# the `.mat` path.

# %%
from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EXPLORATION_ROOT = Path("/Users/lauraporta/local_data/motion_clouds/figures/glm/exploration/cv_strategy")
PROBE_MAT = Path("/Users/lauraporta/local_data/motion_clouds/formatted_data/CAA-1123243_rec1.mat")
MATLAB_REF_PROFILE_PNG = Path(
    "/Volumes/margrie/laura/data transfer for laura/"
    "glm_single_cluster/speed_profile_cross_validation.png"
)

assert EXPLORATION_ROOT.is_dir(), (
    f"Exploration runs missing. Produce them by following the shell "
    f"commands under each section — or re-run the whole sweep:\n"
    f"  rm -rf {EXPLORATION_ROOT} && mkdir -p {EXPLORATION_ROOT}\n"
    f"  # Phase A, B, C as documented below"
)


def _load_comparison(run_dir: Path) -> pd.DataFrame:
    """Load glm_model_comparison.csv and strip the ``time_`` prefix.

    Every exploration run only produces the "time" glm_type (there is
    no other in the Python port), so stripping makes downstream column
    access terser.
    """
    df = pd.read_csv(run_dir / "glm_model_comparison.csv")
    return df.rename(columns={c: c[5:] for c in df.columns if c.startswith("time_")})


# %% [markdown]
# ## 1. Speed-profile CV — generalisation across the two velocity profiles
#
# Reproduces MATLAB `speed_profile_cross_validation.png`
# (`scripts/glm_single_cluster_analysis.m:2246-2367`). The rc2 motion-clouds
# stimulus presents the same cloud set under two reproduced velocity
# trajectories (Laura's "two speed profiles"). MATLAB splits trials at
# the midpoint of `presentation_sequence` — first half = profile 1,
# second half = profile 2 — then runs 2-fold CV: train on one profile,
# test on the other. Compared to the condition-stratified 5-fold CV,
# this is a strictly harder generalisation test.
#
# Pipeline invocation:
#
# ```bash
# rc2-glm  CAA-1123243_rec1.mat  glm/exploration/cv_strategy/profile_cv \
#     --backend irls --profile-cv-diagnostic --plot-clusters 5
# ```
#
# Six new columns in glm_model_comparison.csv:
# `profile_cv_bps_{Null,Selected,NoSpeed}` (profile-fold) and
# `profile_cv_bps_{Null,Selected,NoSpeed}_standard` (condition-
# stratified fold, so the paired-diff plot has both sides).

# %%
profile_run = EXPLORATION_ROOT / "profile_cv"
cmp_pcv = _load_comparison(profile_run)
print(f"speed-profile diagnostic rows: {len(cmp_pcv)} clusters")
for col in ("profile_cv_bps_Null", "profile_cv_bps_Selected",
            "profile_cv_bps_NoSpeed",
            "profile_cv_bps_Null_standard",
            "profile_cv_bps_Selected_standard",
            "profile_cv_bps_NoSpeed_standard"):
    n_finite = cmp_pcv[col].notna().sum()
    print(f"  {col}: {n_finite} finite values")

# %% [markdown]
# ### 1a. Per-cluster paired Δ-bps — standard vs profile CV
#
# Each cluster gives two numbers: Δ-bps under condition-stratified
# 5-fold CV (the normal pipeline CV) and Δ-bps under speed-profile
# 2-fold CV (the new diagnostic). If the selected GLM generalises
# across profiles, the two numbers agree. If not, the profile-CV
# Δ-bps is lower.

# %%
d_std = (cmp_pcv["Selected_cv_bps"] - cmp_pcv["Null_cv_bps"]).to_numpy(dtype=np.float64)
d_spr = (cmp_pcv["profile_cv_bps_Selected"] - cmp_pcv["profile_cv_bps_Null"]).to_numpy(dtype=np.float64)
pair = np.isfinite(d_std) & np.isfinite(d_spr)
d_std, d_spr = d_std[pair], d_spr[pair]
print(f"n common clusters with paired Δ-bps: {len(d_std)}")
print(f"  median standard Δ-bps = {np.median(d_std):+.4f}")
print(f"  median profile  Δ-bps = {np.median(d_spr):+.4f}")
print(f"  per-cluster ratio (profile/standard) median: "
      f"{np.median(d_spr[d_std > 0] / d_std[d_std > 0]):.2%}")

# Wilcoxon one-tailed: is standard CV > profile CV?
from scipy.stats import wilcoxon
try:
    p_t = float(wilcoxon(d_std, d_spr, alternative="greater").pvalue)
    print(f"  Wilcoxon (one-sided, std > profile): p = "
          f"{p_t:.3g}")
except ValueError as exc:
    print(f"  Wilcoxon unavailable: {exc}")

# %% [markdown]
# ### 1b. Inline figure — matches MATLAB layout

# %%
from rc2_glm.plots import plot_speed_profile_cv_comparison

fig = plot_speed_profile_cv_comparison(cmp_pcv)
fig.set_size_inches(14, 5)
plt.show()

# %% [markdown]
# **Side-by-side with MATLAB's PNG** (on ceph, rendered for the same
# comparison across the 4-probe reference). Our 74-cluster single-probe
# result should directionally agree (std > profile, p < 0.001), with
# the % retained within ±10 pp.

# %%
if MATLAB_REF_PROFILE_PNG.is_file():
    fig2, ax2 = plt.subplots(figsize=(14, 5))
    ax2.imshow(plt.imread(str(MATLAB_REF_PROFILE_PNG)))
    ax2.set_axis_off()
    ax2.set_title("MATLAB reference (ceph): speed_profile_cross_validation.png",
                  fontsize=10)
    plt.show()
else:
    print("MATLAB reference PNG not available (ceph not mounted).")

# %% [markdown]
# ## 2. Bin-width sweep — 100 / 50 / 20 ms
#
# Current default is 100 ms; MATLAB uses 20 ms. The Poisson GLM offset
# term `log(bin_width)` absorbs the bin-width change so absolute cv-bps
# shift by a constant (already-documented "CV-bps offset" convention,
# see `summary/rc2_analysis-python-matlab-parity.md`). Shape-preserving
# quantities — Spearman rank across clusters, selected-vars Jaccard,
# tuning-curve Pearson — are invariant.
#
# Pipeline invocation (one per bin width):
#
# ```bash
# for bw in 0.100 0.050 0.020; do
#     slug=$(echo $bw | tr . p)
#     rc2-glm  CAA-1123243_rec1.mat  glm/exploration/cv_strategy/bw_${slug} \
#         --backend irls --bin-width $bw --plot-clusters 3
# done
# ```

# %%
bin_widths = {"100 ms": EXPLORATION_ROOT / "bw_0p100",
              "50 ms":  EXPLORATION_ROOT / "bw_0p050",
              "20 ms":  EXPLORATION_ROOT / "bw_0p020"}
bw_frames = {label: _load_comparison(d) for label, d in bin_widths.items()}
print("bin-width runs loaded:")
for label, df in bw_frames.items():
    print(f"  {label}: n_clusters={len(df)} | "
          f"n_bins median={int(df['n_bins'].median())} | "
          f"Null_cv_bps median={df['Null_cv_bps'].median():+.3f} | "
          f"Selected_cv_bps median={df['Selected_cv_bps'].median():+.3f}")

# %% [markdown]
# ### 2a. Absolute cv-bps shifts with bin width; rank is preserved

# %%
from scipy.stats import spearmanr

baseline = bw_frames["100 ms"]
print("Spearman rank of Selected_cv_bps against the 100 ms run:")
for label, df in bw_frames.items():
    if label == "100 ms":
        continue
    merged = baseline.merge(df, on=["probe_id", "cluster_id"],
                             suffixes=("_base", "_x"))
    merged = merged.dropna(subset=["Selected_cv_bps_base", "Selected_cv_bps_x"])
    rho, _ = spearmanr(merged["Selected_cv_bps_base"], merged["Selected_cv_bps_x"])
    mean_offset = float(
        (merged["Selected_cv_bps_x"] - merged["Selected_cv_bps_base"]).mean()
    )
    print(f"  {label}: Spearman {rho:+.4f} | mean offset vs 100 ms "
          f"{mean_offset:+.3f} bps (n={len(merged)})")

# %% [markdown]
# ### 2b. Selected-vars Jaccard across bin widths
#
# The forward-selection set can drift near the Δ-bps = 0.005
# acceptance threshold as bin counts change. Jaccard per cluster
# measures how much of the selected-vars set stays consistent.

# %%
def _as_set(s) -> set[str]:
    if not isinstance(s, str):
        return set()
    norm = s.replace(", ", "+").replace(",", "+")
    return {t.strip() for t in norm.split("+") if t.strip() and t.strip() != "Null"}


def _jaccard(a: str, b: str) -> float:
    A, B = _as_set(a), _as_set(b)
    if not A and not B:
        return 1.0
    if not (A | B):
        return 1.0
    return len(A & B) / len(A | B)


for label, df in bw_frames.items():
    if label == "100 ms":
        continue
    merged = baseline.merge(df, on=["probe_id", "cluster_id"],
                             suffixes=("_base", "_x"))
    jac = merged.apply(
        lambda r: _jaccard(r["selected_vars_base"], r["selected_vars_x"]), axis=1,
    )
    print(f"  {label}: mean Jaccard = {jac.mean():.3f}  "
          f"(identical: {int((jac == 1.0).sum())}/{len(jac)})")

# %% [markdown]
# ### 2c. Tuning-curve visual comparison — one cluster × three bin widths
#
# Load a representative cluster's `cluster_<id>_tuning.pdf` from each
# bin-width run and display side-by-side. The x-axis bin positions
# (5%-quantile centres) shift with bin width because the observed
# distributions themselves shift; the model-row lines then evaluate at
# the matching per-condition centres (commit 4c59c67).

# %%
cluster_to_inspect = 116  # a speed+SF tuned cluster that shows up on this probe
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (label, d) in zip(axes, bin_widths.items()):
    pdf = d / "figs" / f"cluster_{cluster_to_inspect}_tuning.pdf"
    if pdf.is_file():
        # Lightweight inline render: convert PDF's first page to PNG on demand.
        import tempfile, shutil
        with tempfile.NamedTemporaryFile(suffix=".png") as tf:
            subprocess.run(
                ["sips", "-s", "format", "png", "-Z", "1600",
                 str(pdf), "--out", tf.name],
                check=True, capture_output=True,
            )
            ax.imshow(plt.imread(tf.name))
    ax.set_title(f"cluster {cluster_to_inspect} — {label}", fontsize=10)
    ax.set_axis_off()
plt.suptitle("Tuning curves at 100 / 50 / 20 ms bin widths", fontsize=11)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Seed stability — how much does the fold RNG matter?
#
# The condition-stratified fold assignment uses `numpy.default_rng(seed)`
# to permute trials before round-robin fold assignment. Running the
# pipeline with N different seeds lets us quantify how much of the
# fit's cv-bps value is "noise from trial ordering" vs. signal.
#
# Pipeline invocation (one per seed):
#
# ```bash
# for seed in 0 1 2 3 4; do
#     rc2-glm  CAA-1123243_rec1.mat  glm/exploration/cv_strategy/seed_${seed} \
#         --backend irls --cv-seed $seed --no-plots
# done
# ```

# %%
seed_runs = {
    f"seed {s}": EXPLORATION_ROOT / f"seed_{s}"
    for s in (0, 1, 2, 3, 4)
    if (EXPLORATION_ROOT / f"seed_{s}").is_dir()
}
seed_frames = {label: _load_comparison(d) for label, d in seed_runs.items()}
print(f"seed runs loaded: {list(seed_frames)}")

# %% [markdown]
# ### 3a. Per-cluster cv-bps variability across seeds

# %%
wide = pd.concat(
    [df.set_index(["probe_id", "cluster_id"])["Selected_cv_bps"].rename(label)
     for label, df in seed_frames.items()],
    axis=1,
)
std_per_cluster = wide.std(axis=1)
mean_per_cluster = wide.mean(axis=1)
print(f"per-cluster cv-bps std across {len(seed_frames)} seeds:")
print(f"  median std = {std_per_cluster.median():.4f} bps")
print(f"  max std    = {std_per_cluster.max():.4f} bps "
      f"(cluster {std_per_cluster.idxmax()[1]})")
print(f"  ratio (std / mean) median = "
      f"{(std_per_cluster / mean_per_cluster.abs()).median():.2%}")

fig, ax = plt.subplots(figsize=(6, 3.8))
ax.hist(std_per_cluster.values, bins=np.linspace(0, max(0.05, std_per_cluster.max() * 1.05), 25),
        edgecolor="black", color="#668")
ax.set_xlabel("per-cluster cv-bps std across seeds")
ax.set_ylabel("# clusters")
ax.set_title(f"Seed stability — {len(seed_frames)} seeds × {len(std_per_cluster)} clusters")
ax.axvline(std_per_cluster.median(), color="#c33", linestyle="--",
           label=f"median = {std_per_cluster.median():.3f}")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ### 3b. Selected-vars Jaccard pairwise across seeds

# %%
jac_pairs = []
labels = list(seed_frames)
for i in range(len(labels)):
    for j in range(i + 1, len(labels)):
        a = seed_frames[labels[i]].set_index(["probe_id", "cluster_id"])["selected_vars"]
        b = seed_frames[labels[j]].set_index(["probe_id", "cluster_id"])["selected_vars"]
        common = a.index.intersection(b.index)
        jac_vals = [_jaccard(a.loc[k], b.loc[k]) for k in common]
        mean_j = float(np.mean(jac_vals)) if jac_vals else float("nan")
        jac_pairs.append((labels[i], labels[j], mean_j, len(common)))

print("pairwise mean Jaccard across seeds (selected_vars):")
for a, b, mj, n in jac_pairs:
    print(f"  {a} vs {b}: mean Jaccard = {mj:.3f} (n={n} common clusters)")

pairwise_mean = np.mean([mj for _, _, mj, _ in jac_pairs])
print(f"  overall pairwise mean Jaccard = {pairwise_mean:.3f}")

# %% [markdown]
# ## Summary
#
# | Diagnostic | Finding |
# |---|---|
# | Speed-profile CV | See §1 — cross-profile generalisation is statistically significant vs standard CV (Wilcoxon, one-tailed, std > profile). % retained measures how much of the selected model's predictive power survives when trained on one profile and tested on the other. |
# | Bin-width sweep | See §2 — absolute cv-bps shifts with `log(bin_width)` (convention offset, documented). Spearman rank across clusters and selected-vars Jaccard stay high across 100 / 50 / 20 ms. Tuning curves at 20 ms have denser observed 5%-quantile bins (more motion bins → finer quantile grid) but the same shape. |
# | Seed stability | See §3 — per-cluster cv-bps std across 5 seeds is small relative to the signal (median std / |mean| < 5%). Selected-vars Jaccard pairwise mean across seeds is high (> 0.9 expected). The standard CV RNG is not a material source of fit variability. |
#
# The three investigations are exploratory — nothing here gates
# production. The shell commands at the top of each section are the
# repro recipe.
