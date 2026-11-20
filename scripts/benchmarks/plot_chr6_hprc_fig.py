#!/usr/bin/env python3
"""Generate BioVis-compliant Fig 6: Horizon vs HPRC Minigraph-Cactus on GRCh38.chr6 x CHM13.chr6."""
import json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# BioVis typography & rendering defaults (from richard-burhans/writing/papers/scripts/figure_template.py)
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#6B7280",
    "light_gray": "#D1D5DB",
}

WIDTH_2COL = 170.0 * 0.0393701  # 6.69 inches

OUT_DIR = Path(__file__).resolve().parents[2] / "manuscript" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

data = np.load("/tmp/chr6_comparison_data.npz")
with open("/tmp/chr6_comparison_summary.json") as f:
    summary = json.load(f)

x_hprc = data["x_hprc_bp"] / 1e6
x_hor = data["x_hor_bp"] / 1e6
x_exact = data["x_exact_bp"] / 1e6
truth_segs = data["truth_segs"]
hor_blocks = data["hor_blocks"]

fig = plt.figure(figsize=(WIDTH_2COL, 4.7))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.15, 1.0], hspace=0.40, wspace=0.30)

def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# --- Panel A: Whole-chromosome synteny dot-plot ---
ax_a = fig.add_subplot(gs[0, 0])
# Subsample/plot truth segments >= 500 bp as light gray/sky background, Horizon blocks >= 500 bp in blue
tx_mid = (truth_segs[:, 0] + truth_segs[:, 1]) / 2.0 / 1e6
ty_mid = ((truth_segs[:, 2] + truth_segs[:, 3]) / 2.0 + summary["y_offset_bp"]) / 1e6
ax_a.scatter(tx_mid, ty_mid, s=1.2, c=OI["sky"], alpha=0.35, edgecolors="none", rasterized=True, label="HPRC 1-to-1 segs (≥500 bp)")

hx_mid = (hor_blocks[:, 0] + hor_blocks[:, 1]) / 2.0 / 1e6
hy_mid = ((hor_blocks[:, 2] + hor_blocks[:, 3]) / 2.0 + summary["y_offset_bp"]) / 1e6
ax_a.scatter(hx_mid, hy_mid, s=0.9, c=OI["blue"], alpha=0.55, edgecolors="none", rasterized=True, label="Horizon blocks (≥500 bp)")

ax_a.axvline(61.0, color=OI["orange"], linestyle="--", linewidth=0.9, alpha=0.8)
ax_a.text(63.5, 72.0, "q-arm\nstart", fontsize=6.5, color=OI["vermilion"], va="bottom")

ax_a.set_xlim(0, 171)
ax_a.set_ylim(61, 173)
ax_a.set_xlabel("GRCh38.chr6 position (Mbp)", fontsize=8.5)
ax_a.set_ylabel("CHM13.chr6 position (Mbp)", fontsize=8.5)
ax_a.tick_params(labelsize=7.5)
ax_a.legend(loc="upper left", frameon=False, fontsize=6.5, markerscale=3.0, handletextpad=0.3)
despine(ax_a)
ax_a.text(-0.18, 1.08, "A", transform=ax_a.transAxes, fontsize=11, fontweight="bold", va="top")

# --- Panel B: Coordinate precision tolerances ---
ax_b = fig.add_subplot(gs[0, 1])
labels = ["Exact\n(0 bp)", "≤ 10 bp", "≤ 100 bp", "≤ 1 kb"]
vals = [
    100.0 * summary["exact_coord_precision_on_shared"],
    100.0 * summary["within_10bp_pairs"] / summary["pairs_on_hprc_1to1_x"],
    100.0 * summary["within_100bp_precision_on_shared"],
    100.0 * summary["within_1kb_precision_on_shared"],
]
y_pos = np.arange(len(labels))
bars = ax_b.barh(y_pos, vals, height=0.52, color=[OI["blue"], OI["sky"], OI["green"], OI["gray"]], alpha=0.9)
ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(labels, fontsize=7.5)
ax_b.set_ylim(-0.55, 3.65)
ax_b.invert_yaxis()
ax_b.set_xlim(99.0, 100.12)
ax_b.set_xlabel("Coordinate precision vs HPRC graph (%)", fontsize=8.5)
ax_b.tick_params(labelsize=7.5)
for bar, v in zip(bars, vals):
    ax_b.text(v - 0.03, bar.get_y() + bar.get_height() / 2, f"{v:.2f}%",
              ha="right", va="center", fontsize=7.0, fontweight="bold", color="white")
ax_b.set_title(
    f"n = {summary['pairs_on_hprc_1to1_x']/1e6:.2f}M shared pairs  |  Exact recall: {100*summary['exact_coord_recall_on_hprc_1to1']:.1f}%",
    fontsize=7.0, fontweight="normal", color=OI["black"], loc="left", pad=4,
)
despine(ax_b)
ax_b.text(-0.20, 1.08, "B", transform=ax_b.transAxes, fontsize=11, fontweight="bold", va="top")

# --- Panel C: 1-Mbp binned coverage & exact concordance across GRCh38.chr6 ---
ax_c = fig.add_subplot(gs[1, :])
bins_mb = np.arange(len(x_hprc)) + 0.5
ax_c.fill_between(bins_mb, 0, x_hprc, color=OI["light_gray"], alpha=0.55, label="HPRC 1-to-1 shared bp")
ax_c.plot(bins_mb, x_hor, color=OI["sky"], linewidth=1.1, linestyle="-", label="Horizon MAF aligned bp")
ax_c.plot(bins_mb, x_exact, color=OI["blue"], linewidth=1.2, linestyle="-", label="Exact match (Δy = 0 bp)")

ax_c.axvspan(58.0, 62.0, color=OI["orange"], alpha=0.14, label="Centromere / q-arm")
ax_c.set_xlim(0, 171)
ax_c.set_ylim(0, 1.05)
ax_c.set_xlabel("GRCh38.chr6 position (1 Mbp bins, n = 171)", fontsize=8.5)
ax_c.set_ylabel("Aligned fraction per Mbp", fontsize=8.5)
ax_c.tick_params(labelsize=7.5)
ax_c.legend(loc="upper left", bbox_to_anchor=(0.01, 0.96), ncol=1, frameon=False, fontsize=6.6, handlelength=1.5)
despine(ax_c)
ax_c.text(-0.08, 1.08, "C", transform=ax_c.transAxes, fontsize=11, fontweight="bold", va="top")

for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT_DIR / f"fig6_chr6_hprc_comparison.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("Saved fig6_chr6_hprc_comparison.{png,svg,pdf} to", OUT_DIR)
