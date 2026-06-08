"""Cross-condition transition plots for the split-by-condition GLMs.

Reads the canonical 10-seed run
(``figures/glm/current_splitted_by_condition/``) and renders, per the
2026-06-08 request, plots that show how each cluster's forward-selected
best model changes between V, T_Vstatic and VT.

Distinction worth keeping straight (these plots use the FIRST sense):
  * "Additive" model = the rc2_glm GLM with all main effects summed and
    no interaction terms (Hardcastle sense). Used here only as the source
    of a complete per-variable marginal tuning curve for every cluster.
  * NOT the V+T=VT complete-additive decomposition (that is the
    trial-averaged speed_two_models analysis — a different question).

Figures (into ``.../figs/transitions/``):
  1. selection_fraction_by_condition.pdf
       grouped bars: fraction of the 88 selecting each variable, per condition.
  2. speed_tf_transition_{V_to_VT,T_to_VT}.pdf
       4×4 transition heatmaps over the {none, TF, Speed, Speed+TF}
       Speed/TF selection category — the reattribution at a glance.
  3. tf_amplitude_V_vs_VT.pdf / speed_amplitude_T_vs_VT.pdf
       paired scatters of marginal-tuning amplitude (max-min log-gain of
       the Additive model's curve), y=x diagonal, points coloured by
       whether the cluster swapped TF→Speed (V→VT).
  4. cvbps_delta_pairs.pdf
       Δ cv_bps (Selected-Null) paired V-vs-VT and T-vs-VT.

Usage::  python -m scripts.plot_split_by_condition_transitions
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("split_transitions")

HOME = Path.home()
RUN = HOME / "local_data" / "motion_clouds" / "figures" / "glm" / "current_splitted_by_condition"
FIGDIR = RUN / "figs" / "transitions"
PROBES = ("CAA-1123243_rec1", "CAA-1123244_rec1", "CAA-1123466_rec1", "CAA-1123467_rec1")
CONDITIONS = ("V", "T_Vstatic", "VT")
KEY = ["probe_id", "cluster_id"]


def load_comparison(cond: str) -> pd.DataFrame:
    return pd.read_csv(RUN / f"glm_model_comparison_{cond}.csv")


def speed_tf_category(df: pd.DataFrame) -> pd.Series:
    """Compact Speed/TF selection label per cluster: none/TF/Speed/Speed+TF."""
    sp = df.get("time_is_speed_tuned", 0).astype(bool)
    tf = df.get("time_is_tf_tuned", 0).astype(bool)
    out = np.where(sp & tf, "Speed+TF",
          np.where(sp & ~tf, "Speed",
          np.where(~sp & tf, "TF", "none")))
    return pd.Series(out, index=df.index)


# ----------------------------------------------------------------------
# 1. selection fraction per condition
# ----------------------------------------------------------------------
def plot_selection_fraction() -> None:
    variables = [("Speed", "time_is_speed_tuned"), ("TF", "time_is_tf_tuned"),
                 ("SF", "time_is_sf_tuned"), ("OR", "time_is_or_tuned")]
    data = {}
    n = None
    for cond in CONDITIONS:
        df = load_comparison(cond)
        n = len(df)
        data[cond] = [100.0 * df[col].sum() / n if col in df else 0.0
                      for _, col in variables]
    x = np.arange(len(variables))
    w = 0.26
    fig, ax = plt.subplots(figsize=(7, 4))
    colours = {"V": "tab:blue", "T_Vstatic": "tab:green", "VT": "tab:red"}
    for i, cond in enumerate(CONDITIONS):
        bars = ax.bar(x + (i - 1) * w, data[cond], w, label=cond, color=colours[cond])
        for b, v in zip(bars, data[cond]):
            if v > 0:
                ax.text(b.get_x() + b.get_width() / 2, v + 0.8, f"{v:.0f}",
                        ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([v for v, _ in variables])
    ax.set_ylabel(f"% of cohort selecting (n={n})")
    ax.set_title("Forward-selected variables per condition (10-seed 7/10)")
    ax.legend(title="condition")
    _save(fig, "selection_fraction_by_condition")


# ----------------------------------------------------------------------
# 2. Speed/TF transition heatmaps
# ----------------------------------------------------------------------
def plot_transition(cond_from: str, cond_to: str) -> None:
    a = load_comparison(cond_from)[KEY].assign(cat_from=speed_tf_category(load_comparison(cond_from)))
    b = load_comparison(cond_to)[KEY].assign(cat_to=speed_tf_category(load_comparison(cond_to)))
    m = a.merge(b, on=KEY)
    order = ["none", "TF", "Speed", "Speed+TF"]
    tab = pd.crosstab(m["cat_from"], m["cat_to"]).reindex(index=order, columns=order, fill_value=0)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(tab.values, cmap="Blues", aspect="equal")
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    ax.set_xlabel(f"{cond_to} (to)"); ax.set_ylabel(f"{cond_from} (from)")
    for i in range(len(order)):
        for j in range(len(order)):
            v = int(tab.values[i, j])
            if v:
                ax.text(j, i, str(v), ha="center", va="center",
                        color="white" if v > tab.values.max() / 2 else "black", fontsize=10)
    ax.set_title(f"Speed/TF selection transition\n{cond_from} → {cond_to}")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="clusters")
    _save(fig, f"speed_tf_transition_{cond_from}_to_{cond_to}")


# ----------------------------------------------------------------------
# 3. tuning-amplitude pair plots
# ----------------------------------------------------------------------
def tuning_amplitude(cond: str, variable: str) -> pd.DataFrame:
    """max-min log-gain of the Additive-model marginal curve, per cluster."""
    rows = []
    for probe in PROBES:
        csv = RUN / "_runs" / probe / cond / "diagnostics" / "tuning_curves.csv"
        if not csv.exists():
            continue
        t = pd.read_csv(csv)
        sub = t[(t["model"] == "Additive") & (t["variable"] == variable)]
        for cid, g in sub.groupby("cluster_id"):
            rows.append({"probe_id": probe, "cluster_id": int(cid),
                         "amp": float(g["log_gain"].max() - g["log_gain"].min())})
    return pd.DataFrame(rows)


def plot_amplitude_pair(variable: str, cond_a: str, cond_b: str, swap_mask_fn) -> None:
    A = tuning_amplitude(cond_a, variable).rename(columns={"amp": "amp_a"})
    B = tuning_amplitude(cond_b, variable).rename(columns={"amp": "amp_b"})
    m = A.merge(B, on=KEY)
    swap = swap_mask_fn(m)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(m.loc[~swap, "amp_a"], m.loc[~swap, "amp_b"], s=22,
               c="0.6", label="other", edgecolor="none")
    ax.scatter(m.loc[swap, "amp_a"], m.loc[swap, "amp_b"], s=34,
               c="tab:red", label="TF→Speed swap (V→VT)", edgecolor="k", linewidth=0.4)
    hi = float(np.nanmax([m["amp_a"].max(), m["amp_b"].max()])) * 1.05
    ax.plot([0, hi], [0, hi], ls="--", c="0.5", lw=1)
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    ax.set_xlabel(f"{variable} tuning amplitude — {cond_a}  (Δ log-gain)")
    ax.set_ylabel(f"{variable} tuning amplitude — {cond_b}  (Δ log-gain)")
    ax.set_title(f"{variable} marginal-tuning amplitude: {cond_a} vs {cond_b}\n(Additive model)")
    ax.legend(fontsize=8)
    _save(fig, f"{variable.lower()}_amplitude_{cond_a}_vs_{cond_b}")


# ----------------------------------------------------------------------
# 4. cv_bps delta pairs
# ----------------------------------------------------------------------
def plot_cvbps_pairs() -> None:
    col = "time_delta_selected_vs_null"
    dfs = {c: load_comparison(c)[KEY + [col]].rename(columns={col: c}) for c in CONDITIONS}
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, (ca, cb) in zip(axes, [("V", "VT"), ("T_Vstatic", "VT")]):
        m = dfs[ca].merge(dfs[cb], on=KEY)
        ax.scatter(m[ca], m[cb], s=24, c="tab:purple", edgecolor="none", alpha=0.8)
        hi = float(np.nanmax([m[ca].max(), m[cb].max()])) * 1.05
        lo = min(0.0, float(np.nanmin([m[ca].min(), m[cb].min()])))
        ax.plot([lo, hi], [lo, hi], ls="--", c="0.5", lw=1)
        ax.axhline(0, c="0.8", lw=0.8); ax.axvline(0, c="0.8", lw=0.8)
        ax.set_xlabel(f"Δ cv_bps (Sel−Null) — {ca}")
        ax.set_ylabel(f"Δ cv_bps (Sel−Null) — {cb}")
        ax.set_title(f"{ca} vs {cb}")
    fig.suptitle("Selected-model gain over baseline, paired per cluster", fontweight="bold")
    _save(fig, "cvbps_delta_pairs")


# ----------------------------------------------------------------------
def _swap_tf_to_speed(merged: pd.DataFrame) -> pd.Series:
    """clusters TF-selected in V that lose TF and gain Speed in VT."""
    V = load_comparison("V")[KEY].assign(tfV=speed_tf_category(load_comparison("V")).isin(["TF", "Speed+TF"]))
    VT = load_comparison("VT")
    VTc = VT[KEY].assign(tfVT=VT["time_is_tf_tuned"].astype(bool),
                         spVT=VT["time_is_speed_tuned"].astype(bool))
    j = V.merge(VTc, on=KEY)
    j["swap"] = j["tfV"] & (~j["tfVT"]) & j["spVT"]
    out = merged.merge(j[KEY + ["swap"]], on=KEY, how="left")["swap"].fillna(False)
    out.index = merged.index
    return out


def plot_condition_alluvial() -> None:
    """Alluvial of the Speed/TF selected category across T_Vstatic → V → VT.

    Each cluster is a ribbon flowing through the three conditions; ribbon
    colour = its VT (integration) category. Shows how the forward-selected
    model transitions across conditions (e.g. Speed-in-T + TF-in-V → what in VT).
    """
    cats = ["none", "TF", "Speed", "Speed+TF"]
    cidx = {c: i for i, c in enumerate(cats)}
    pal = {"none": "0.8", "TF": "tab:orange", "Speed": "tab:green",
           "Speed+TF": "tab:purple"}
    cols = ["T_Vstatic", "V", "VT"]
    catcol = {"T_Vstatic": "catT", "V": "catV", "VT": "catVT"}
    frames = {}
    for cond in cols:
        d = load_comparison(cond)
        frames[cond] = pd.DataFrame({**{k: d[k] for k in KEY},
                                     catcol[cond]: speed_tf_category(d).values})
    m = (frames["T_Vstatic"].merge(frames["V"], on=KEY)
         .merge(frames["VT"], on=KEY))
    groups = m.groupby(["catT", "catV", "catVT"]).size().reset_index(name="n")
    xpos = {"T_Vstatic": 0.0, "V": 1.0, "VT": 2.0}
    GAP, BOXW = 3.0, 0.10

    # y-offset of each group within each column (category boxes contiguous).
    ypos: dict = {}
    box_extent: dict = {}  # (col,cat) -> [ymin,ymax]
    for col in cols:
        cc = catcol[col]
        order = sorted(range(len(groups)), key=lambda i: (
            cidx[groups.loc[i, cc]], cidx[groups.loc[i, "catVT"]],
            cidx[groups.loc[i, "catT"]], cidx[groups.loc[i, "catV"]]))
        y, prev = 0.0, None
        for i in order:
            c = groups.loc[i, cc]
            if prev is not None and c != prev:
                y += GAP
            ypos[(col, i)] = y
            ext = box_extent.setdefault((col, c), [y, y])
            ext[0] = min(ext[0], y); ext[1] = max(ext[1], y + groups.loc[i, "n"])
            y += groups.loc[i, "n"]
            prev = c

    fig, ax = plt.subplots(figsize=(8.5, 6))
    # Ribbons first (grey 'none' underneath, coloured on top), then boxes.
    draw_order = sorted(range(len(groups)),
                        key=lambda i: cidx[groups.loc[i, "catVT"]])
    for i in draw_order:
        n = int(groups.loc[i, "n"]); color = pal[groups.loc[i, "catVT"]]
        alpha = 0.25 if groups.loc[i, "catVT"] == "none" else 0.7
        for a, b in (("T_Vstatic", "V"), ("V", "VT")):
            xa, xb = xpos[a] + BOXW, xpos[b]
            ya, yb = ypos[(a, i)], ypos[(b, i)]
            ax.add_patch(plt.Polygon([(xa, ya), (xb, yb), (xb, yb + n), (xa, ya + n)],
                                     closed=True, color=color, alpha=alpha, lw=0))
    for i in range(len(groups)):
        n = int(groups.loc[i, "n"]); color = pal[groups.loc[i, "catVT"]]
        for col in cols:
            ax.add_patch(plt.Rectangle((xpos[col], ypos[(col, i)]), BOXW, n,
                                       color=color, alpha=1.0, lw=0))
    for (col, c), (lo, hi) in box_extent.items():
        ax.text(xpos[col] + BOXW / 2, (lo + hi) / 2, f"{c}\n({int(hi - lo)})",
                ha="center", va="center", fontsize=7, fontweight="bold")
    ax.set_xticks(list(xpos.values()))
    ax.set_xticklabels(["T_Vstatic", "V", "VT"], fontsize=11)
    ax.set_yticks([]); ax.set_ylabel("clusters")
    ax.set_xlim(-0.3, 2.4)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=pal[c]) for c in cats]
    ax.legend(handles, [f"VT: {c}" for c in cats], fontsize=8, loc="upper right",
              title="ribbon colour = VT category", title_fontsize=8)
    ax.set_title("Selected Speed/TF model across conditions  (T_Vstatic → V → VT)\n"
                 "ribbons = clusters, coloured by their VT integration category",
                 fontsize=11, fontweight="bold")
    _save(fig, "condition_alluvial_T_V_VT")


def plot_condition_integration() -> None:
    """How T (Speed?) × V (TF?) selections combine into the VT selected model.

    T_Vstatic can only select Speed; V can only select TF. Group clusters by
    (Speed-in-T, TF-in-V) — the four input combinations — and stack the VT
    category within each, to show the T→V→VT transition / integration cleanly
    (no alluvial muddiness). E.g. cells Speed-in-T AND TF-in-V → do they become
    Speed+TF in VT, or does one modality win?
    """
    frames = {}
    cc = {"T_Vstatic": "catT", "V": "catV", "VT": "catVT"}
    for cond in ("T_Vstatic", "V", "VT"):
        d = load_comparison(cond)
        frames[cond] = pd.DataFrame({**{k: d[k] for k in KEY},
                                     cc[cond]: speed_tf_category(d).values})
    m = frames["T_Vstatic"].merge(frames["V"], on=KEY).merge(frames["VT"], on=KEY)
    m["T_speed"] = m["catT"].isin(["Speed", "Speed+TF"])
    m["V_tf"] = m["catV"].isin(["TF", "Speed+TF"])
    combos = [(False, False), (False, True), (True, False), (True, True)]
    labels = ["T:–\nV:–", "T:–\nV:TF", "T:Speed\nV:–", "T:Speed\nV:TF"]
    vt_cats = ["none", "TF", "Speed", "Speed+TF"]
    pal = {"none": "0.8", "TF": "tab:orange", "Speed": "tab:green",
           "Speed+TF": "tab:purple"}
    fig, ax = plt.subplots(figsize=(8, 5))
    for gi, (ts, vt) in enumerate(combos):
        sub = m[(m["T_speed"] == ts) & (m["V_tf"] == vt)]
        bottom = 0
        for c in vt_cats:
            cnt = int((sub["catVT"] == c).sum())
            if cnt:
                ax.bar(gi, cnt, bottom=bottom, color=pal[c], edgecolor="white",
                       width=0.7)
                ax.text(gi, bottom + cnt / 2, str(cnt), ha="center", va="center",
                        fontsize=8, color="white" if c != "none" else "black")
            bottom += cnt
        ax.text(gi, bottom + 0.5, f"n={len(sub)}", ha="center", va="bottom",
                fontsize=8, fontweight="bold")
    ax.set_xticks(range(len(combos)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlabel("input selection (T_Vstatic × V)")
    ax.set_ylabel("clusters")
    handles = [plt.Rectangle((0, 0), 1, 1, color=pal[c]) for c in vt_cats]
    ax.legend(handles, [f"VT: {c}" for c in vt_cats], title="VT selected model",
              fontsize=8, title_fontsize=8)
    ax.set_title("T × V selections → VT integration\n"
                 "(stack = VT Speed/TF category within each input combination)",
                 fontsize=11, fontweight="bold")
    _save(fig, "condition_integration_T_V_to_VT")


def _save(fig, name: str) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    path = FIGDIR / f"{name}.pdf"
    fig.savefig(path)
    fig.savefig(FIGDIR / f"{name}.png", dpi=120)
    plt.close(fig)
    log.info("wrote %s", path)


def main() -> int:
    plot_condition_integration()
    plot_condition_alluvial()
    plot_selection_fraction()
    plot_transition("T_Vstatic", "V")
    plot_transition("V", "VT")
    plot_transition("T_Vstatic", "VT")
    plot_amplitude_pair("TF", "V", "VT", _swap_tf_to_speed)
    plot_amplitude_pair("Speed", "T_Vstatic", "VT", _swap_tf_to_speed)
    plot_cvbps_pairs()
    log.info("done — figures under %s", FIGDIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
