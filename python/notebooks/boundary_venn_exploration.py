"""Δ-bps selection-boundary × activity + 4-probe {Speed,TF,SF,OR} UpSet.

CSV-only exploration of the 4-probe ``glm_out_all/`` run — answers two
questions about the current forward-selection output without re-fitting:

1. **Boundary × activity.** Is the Δ-bps = 0.005 selection threshold
   mainly fragile for low-spike-count clusters? (Headline: no — Q1 is
   the LEAST fragile; Q2/Q3 peak. TF decisions dominate the fragile
   rows; Speed is almost never marginal.)
2. **4-probe UpSet of {Speed, TF, SF, OR}** on the 88 common-filtered
   clusters. Nine non-empty cells; ``OR`` never participates in a 3+
   way intersection; 6 "(none)" clusters concentrate in probes 1123466
   and 1123467.

Inputs (read-only):
  - ``glm_out_all/glm_selection_history.csv``
  - ``glm_out_all/glm_model_comparison.csv``

Outputs (written into ``OUT_DIR``):
  - ``boundary_vs_activity.pdf`` — three-panel: per-round scatter,
    quartile rate bar, candidate-identity breakdown.
  - ``upset_4way.pdf``           — UpSet plot with set-size sidebar.
  - ``venn_per_probe_stacked.pdf`` — per-probe subset composition.
  - CSVs: ``boundary_clusters.csv``, ``fragile_by_candidate.csv``,
    ``venn_cells.csv``, ``venn_cells_by_probe.csv``.

Companion vault note (with interpretation):
``results/motion-clouds/2026-04-24-boundary-venn.md`` in
``lauraporta/obsidian-agents``.

Same jupytext percent format as ``parity_report.py`` and
``cv_strategy_exploration.py`` (``# %%`` separates cells,
``# %% [markdown]`` marks markdown).
"""

# %% [markdown]
# # Boundary × activity, and 4-probe {Speed,TF,SF,OR} UpSet
#
# Written 2026-04-24 after the "full-trust" milestone (20/20 parity
# gates, Jaccard 0.932, tuning-curve Pearson frac=0.855/0.900 on the 4
# passive probes). Answers two follow-up questions posed by Laura while
# looking at `findings-so-far.md` §5 (10% of selections flip on fold
# noise):
#
#   Q1. "Was Δ-bps = 0.005 mainly a big deal for clusters that were
#       not that active?"  → **No.** Q1 (lowest spikes) is the LEAST
#       fragile. Q2/Q3 peak. TF dominates fragile decisions; Speed
#       rarely marginal.
#
#   Q2. "Give me a 4-probe Venn of {Speed, TF, SF, OR} on the filtered
#       clusters." → 9 non-empty cells; OR never in 3+ way
#       intersection.

# %%
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, spearmanr

# Inputs: the current 4-probe run on glm-improvements.
GLM_DIR = Path(
    "/Users/lauraporta/local_data/motion_clouds/figures/glm_out_all"
)
# Outputs: parallel directory alongside glm_out_all/ — does not collide
# with the pipeline's own validation/ subtree.
OUT_DIR = Path(
    "/Users/lauraporta/local_data/motion_clouds/figures/boundary_venn"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Selection threshold used by forward_selection (config default).
THRESH = 0.005
# A row is "fragile" if its delta_bps sits in a ±band straddling 0.005 —
# a CV-fold/seed perturbation could flip accept↔reject. Band chosen
# empirically from the observed current-run distribution:
#   accepted rows  : delta_bps ∈ [0.0052, 0.745]  (23 rows in [0.005, 0.010])
#   rejected rows  : delta_bps ∈ [-2.86, 0.0047]  (38 rows in [0.001, 0.005])
FRAGILE_LOW = 0.002
FRAGILE_HIGH = 0.008


# %%
# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------
hist = pd.read_csv(GLM_DIR / "glm_selection_history.csv")
mcmp = pd.read_csv(GLM_DIR / "glm_model_comparison.csv")

print(
    f"selection_history: {len(hist):>4d} rows, "
    f"{hist[['probe_id','cluster_id']].drop_duplicates().shape[0]} "
    "unique clusters"
)
print(f"model_comparison : {len(mcmp):>4d} rows (one per retained cluster)")

# %% [markdown]
# ## Part 1 — Δ-bps = 0.005 × cluster activity
#
# A row is fragile if ``delta_bps ∈ (0.002, 0.008)``. A cluster is
# fragile if ANY of its selection rounds lands in that band. Tracks
# both "barely accepted" (`added=True`, δ ∈ [0.005, 0.008]) and
# "barely rejected" (`added=False`, δ ∈ [0.002, 0.005]).

# %%
hist["is_fragile_row"] = (
    (hist["delta_bps"] > FRAGILE_LOW) & (hist["delta_bps"] < FRAGILE_HIGH)
)
hist["fragility_kind"] = np.where(
    hist["is_fragile_row"],
    np.where(hist["added"], "fragile_accept", "fragile_reject"),
    "non_fragile",
)
hist["dist_to_thresh"] = (hist["delta_bps"] - THRESH).abs()

per_cluster = (
    hist.sort_values("dist_to_thresh")
    .groupby(["probe_id", "cluster_id"], as_index=False)
    .agg(
        nearest_delta=("delta_bps", "first"),
        nearest_added=("added", "first"),
        nearest_candidate=("best_candidate", "first"),
        nearest_dist=("dist_to_thresh", "first"),
        any_fragile=("is_fragile_row", "any"),
        n_fragile_rows=("is_fragile_row", "sum"),
    )
)

act = (
    mcmp[
        [
            "probe_id",
            "cluster_id",
            "time_n_spikes",
            "time_n_bins",
            "n_trials",
            "time_n_selected_vars",
            "time_selected_vars",
        ]
    ]
    .rename(columns={"n_trials": "time_n_trials"})
    .copy()
)
act["spikes_per_bin"] = act["time_n_spikes"] / act["time_n_bins"]

boundary = per_cluster.merge(act, on=["probe_id", "cluster_id"], how="inner")
boundary["at_boundary"] = boundary["any_fragile"]
print(
    f"boundary-merged  : {len(boundary)} clusters, "
    f"{boundary.at_boundary.sum()} flagged fragile"
)

# %%
boundary["activity_quartile"] = pd.qcut(
    boundary["spikes_per_bin"],
    q=4,
    labels=["Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"],
)
rate_tbl = boundary.groupby("activity_quartile", observed=True).agg(
    n_clusters=("cluster_id", "count"),
    n_at_boundary=("at_boundary", "sum"),
    median_spikes_per_bin=("spikes_per_bin", "median"),
    median_n_spikes=("time_n_spikes", "median"),
)
rate_tbl["pct_at_boundary"] = (
    100 * rate_tbl["n_at_boundary"] / rate_tbl["n_clusters"]
).round(1)
print("\n=== Δ-bps near-boundary fraction by activity quartile ===")
print(rate_tbl.to_string())

rho, pval = spearmanr(boundary["time_n_spikes"], boundary["nearest_dist"])
print(f"\nSpearman ρ  (n_spikes vs dist-to-0.005): {rho:+.3f}  (p={pval:.3g})")
rho2, pval2 = spearmanr(boundary["spikes_per_bin"], boundary["nearest_dist"])
print(f"Spearman ρ  (spikes/bin vs dist-to-0.005): {rho2:+.3f}  (p={pval2:.3g})")

nsp_bdy = boundary.loc[boundary.at_boundary, "time_n_spikes"]
nsp_no = boundary.loc[~boundary.at_boundary, "time_n_spikes"]
u, p_mw = mannwhitneyu(nsp_bdy, nsp_no, alternative="less")
print(
    f"Mann-Whitney one-sided (boundary < non-boundary, n_spikes): "
    f"U={u:.0f}, p={p_mw:.3g}  "
    f"[n_boundary={len(nsp_bdy)}, n_nonboundary={len(nsp_no)}]"
)
print(f"  median n_spikes  boundary    : {nsp_bdy.median():.0f}")
print(f"  median n_spikes  non-boundary: {nsp_no.median():.0f}")

boundary.to_csv(OUT_DIR / "boundary_clusters.csv", index=False)

# %%
fragile_rows = hist[hist.is_fragile_row].copy()
cand_tbl = (
    fragile_rows.groupby(["best_candidate", "added"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=[True, False], fill_value=0)
    .rename(columns={True: "fragile_accept", False: "fragile_reject"})
)
cand_tbl["total"] = cand_tbl.sum(axis=1)
cand_tbl = cand_tbl.sort_values("total", ascending=False)
print("\n=== Fragile rows by candidate variable ===")
print(cand_tbl.to_string())
cand_tbl.to_csv(OUT_DIR / "fragile_by_candidate.csv")

# %%
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

ax = axes[0]
hh = hist.merge(
    act[["probe_id", "cluster_id", "time_n_spikes"]],
    on=["probe_id", "cluster_id"],
)
col_row = np.where(hh["added"], "C0", "C3")
ax.scatter(
    hh["time_n_spikes"], hh["delta_bps"],
    c=col_row, s=18, alpha=0.5, edgecolors="none",
)
ax.axhline(THRESH, color="k", lw=1.0, label=f"Δ-bps = {THRESH}")
ax.axhspan(
    FRAGILE_LOW, FRAGILE_HIGH, color="orange", alpha=0.15,
    label=f"fragile band ({FRAGILE_LOW}–{FRAGILE_HIGH})",
)
ax.set_xscale("log")
ax.set_yscale("symlog", linthresh=1e-3)
ax.set_xlabel("cluster spikes (log)")
ax.set_ylabel("Δ-bps per round (symlog)")
ax.set_title("Per-round Δ-bps vs cluster activity")
ax.legend(
    handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="C0",
               markersize=8, label="accepted (added=True)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="C3",
               markersize=8, label="rejected (added=False)"),
        Line2D([0], [0], color="k", lw=1, label=f"threshold = {THRESH}"),
    ],
    loc="lower right",
    fontsize=8,
)

ax = axes[1]
labels = list(rate_tbl.index)
ax.bar(labels, rate_tbl["pct_at_boundary"].values, color="C2")
for i, v_ in enumerate(rate_tbl["pct_at_boundary"].values):
    ax.text(
        i, v_ + 0.7,
        f"{v_:.1f}%\n(n={rate_tbl['n_clusters'].iloc[i]})",
        ha="center", fontsize=9,
    )
ax.set_ylabel("% fragile clusters")
ax.set_xlabel("Activity quartile (spikes / bin)")
ax.set_title(
    "Fragile-cluster rate by activity quartile\n"
    f"(fragile = Δ-bps ∈ ({FRAGILE_LOW}, {FRAGILE_HIGH}) in any round)"
)
ax.set_ylim(0, max(rate_tbl["pct_at_boundary"].max() * 1.3, 20))
plt.sca(axes[1])
plt.xticks(rotation=15, ha="right")

ax = axes[2]
y = np.arange(len(cand_tbl))
ax.barh(
    y, cand_tbl["fragile_accept"].values,
    color="C0", label="fragile accept (δ∈[0.005,0.008])",
)
ax.barh(
    y, cand_tbl["fragile_reject"].values,
    left=cand_tbl["fragile_accept"].values,
    color="C3", label="fragile reject (δ∈[0.002,0.005])",
)
for yi, (_idx, row) in enumerate(cand_tbl.iterrows()):
    ax.text(row["total"] + 0.2, yi, f"{int(row['total'])}",
            va="center", fontsize=9)
ax.set_yticks(y)
ax.set_yticklabels(cand_tbl.index)
ax.set_xlabel("# fragile rows")
ax.set_title(
    "Fragile decisions by variable identity\n"
    "(TF dominates, Speed rarely marginal)"
)
ax.legend(loc="lower right", fontsize=8, frameon=False)
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(OUT_DIR / "boundary_vs_activity.pdf")
plt.savefig(OUT_DIR / "boundary_vs_activity.png", dpi=150)
print(f"\nwrote {OUT_DIR / 'boundary_vs_activity.pdf'}")

# %% [markdown]
# ## Part 2 — 4-probe UpSet of {Speed, TF, SF, OR}
#
# One row per cluster; boolean membership taken from
# ``time_is_{speed,tf,sf,or}_tuned``.

# %%
v = mcmp[
    [
        "probe_id",
        "cluster_id",
        "time_is_speed_tuned",
        "time_is_tf_tuned",
        "time_is_sf_tuned",
        "time_is_or_tuned",
    ]
].copy()
v.columns = ["probe_id", "cluster_id", "Speed", "TF", "SF", "OR"]

probes = sorted(v["probe_id"].unique())
print("\n=== Variable-selection counts per probe ===")
for p in probes:
    sub = v[v.probe_id == p]
    print(
        f"  {p}: n={len(sub):>3d}   "
        f"Speed={sub.Speed.sum():>3d}   TF={sub.TF.sum():>3d}   "
        f"SF={sub.SF.sum():>3d}   OR={sub.OR.sum():>3d}"
    )
print(
    f"  ALL       : n={len(v):>3d}   "
    f"Speed={v.Speed.sum():>3d}   TF={v.TF.sum():>3d}   "
    f"SF={v.SF.sum():>3d}   OR={v.OR.sum():>3d}"
)


def cell_key(row):
    bits = []
    if row.Speed:
        bits.append("Speed")
    if row.TF:
        bits.append("TF")
    if row.SF:
        bits.append("SF")
    if row.OR:
        bits.append("OR")
    return " ∩ ".join(bits) if bits else "(none)"


v["cell"] = v.apply(cell_key, axis=1)
cell_counts = v["cell"].value_counts().rename("count").to_frame()
cell_counts.index.name = "subset"
print("\n=== 4-way Venn cell counts (across 4 probes) ===")
print(cell_counts.to_string())
cell_counts.to_csv(OUT_DIR / "venn_cells.csv")

per_probe = (
    v.groupby(["probe_id", "cell"])
    .size()
    .unstack(fill_value=0)
    .T.rename_axis("subset", axis=0)
)
per_probe["ALL"] = per_probe.sum(axis=1)
per_probe.to_csv(OUT_DIR / "venn_cells_by_probe.csv")
print("\n=== Per-probe Venn cells ===")
print(per_probe.to_string())

# %%
# UpSet-style plot — clearer than any drawn 4-way Venn.
sets_order = ["Speed", "TF", "SF", "OR"]
cells_sorted = cell_counts.sort_values("count", ascending=False)

fig = plt.figure(figsize=(11, 6.0))
gs = fig.add_gridspec(
    2, 2, width_ratios=[1, 6], height_ratios=[3, 2],
    hspace=0.05, wspace=0.02,
)

ax_bar = fig.add_subplot(gs[0, 1])
x = np.arange(len(cells_sorted))
ax_bar.bar(x, cells_sorted["count"].values, color="#4169E1")
for xi, c in zip(x, cells_sorted["count"].values):
    ax_bar.text(xi, c + 0.4, str(c), ha="center", fontsize=9)
ax_bar.set_xticks(x)
ax_bar.set_xticklabels([])
ax_bar.set_ylabel("# clusters")
ax_bar.set_title(
    f"4-probe UpSet of {{Speed, TF, SF, OR}} — filtered clusters (n={len(v)})"
)
ax_bar.margins(x=0.01)
ax_bar.spines[["top", "right"]].set_visible(False)

ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
for i, name in enumerate(sets_order):
    yrow = len(sets_order) - 1 - i
    for xi, (cell_name, _) in enumerate(cells_sorted.iterrows()):
        filled = name in cell_name
        ax_dot.scatter(
            xi, yrow,
            s=140,
            c="#4169E1" if filled else "#DDDDDD",
            edgecolors="k" if filled else "none",
            linewidths=0.6,
            zorder=3,
        )
for xi, (cell_name, _) in enumerate(cells_sorted.iterrows()):
    members = [name for name in sets_order if name in cell_name]
    if len(members) >= 2:
        ys = [len(sets_order) - 1 - sets_order.index(m) for m in members]
        ax_dot.plot(
            [xi, xi], [min(ys), max(ys)],
            color="#4169E1", lw=2, zorder=2,
        )
ax_dot.set_yticks(range(len(sets_order))[::-1])
ax_dot.set_yticklabels(sets_order)
ax_dot.set_xticks(x)
ax_dot.set_xticklabels([""] * len(x))
ax_dot.set_ylim(-0.5, len(sets_order) - 0.5)
ax_dot.spines[["top", "right", "bottom"]].set_visible(False)
ax_dot.tick_params(axis="x", length=0)
ax_dot.set_xlabel("")

ax_side = fig.add_subplot(gs[1, 0], sharey=ax_dot)
totals = [v[name].sum() for name in sets_order]
ax_side.barh(range(len(sets_order))[::-1], totals, color="#4169E1")
for i, t in enumerate(totals):
    ax_side.text(
        t - 0.5, len(sets_order) - 1 - i, str(t),
        ha="right", va="center", color="white", fontsize=10, weight="bold",
    )
ax_side.invert_xaxis()
ax_side.set_xlabel("Set size")
ax_side.spines[["top", "right"]].set_visible(False)

ax_dot.yaxis.tick_right()
ax_dot.yaxis.set_label_position("right")
for lab in ax_dot.get_yticklabels():
    lab.set_fontsize(11)
    lab.set_weight("bold")

ax_empty = fig.add_subplot(gs[0, 0])
ax_empty.axis("off")

plt.savefig(OUT_DIR / "upset_4way.pdf", bbox_inches="tight")
plt.savefig(OUT_DIR / "upset_4way.png", dpi=150, bbox_inches="tight")
print(f"\nwrote {OUT_DIR / 'upset_4way.pdf'}")

# %%
# Per-probe stacked bar (Venn-as-stack alternative).
fig, ax = plt.subplots(figsize=(8, 5))
per_probe_no_all = per_probe.drop(columns=["ALL"], errors="ignore")
per_probe_no_all = per_probe_no_all.loc[cells_sorted.index]
bottoms = np.zeros(len(per_probe_no_all.columns))
cmap = plt.get_cmap("tab20")
for i, (cell_name, row) in enumerate(per_probe_no_all.iterrows()):
    ax.bar(
        per_probe_no_all.columns, row.values,
        bottom=bottoms, label=cell_name,
        color=cmap(i % 20),
    )
    bottoms = bottoms + row.values
ax.set_ylabel("# clusters")
ax.set_title("Subset composition per probe")
ax.legend(
    bbox_to_anchor=(1.02, 1.0), loc="upper left",
    fontsize=8, frameon=False, ncol=1,
)
plt.tight_layout()
plt.savefig(OUT_DIR / "venn_per_probe_stacked.pdf")
plt.savefig(OUT_DIR / "venn_per_probe_stacked.png", dpi=150)
print(f"wrote {OUT_DIR / 'venn_per_probe_stacked.pdf'}")

print("\nDone.")
