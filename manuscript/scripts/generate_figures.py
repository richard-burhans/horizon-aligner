#!/usr/bin/env python3
"""Generate publication-grade BioVis figures for the Horizon manuscript (Okabe-Ito palette)."""

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colorblind-safe palette
OI_BLUE = "#0072B2"
OI_ORANGE = "#E69F00"
OI_VERMILION = "#D55E00"
OI_GREEN = "#009E73"
OI_SKY = "#56B4E9"
OI_PURPLE = "#CC79A7"
OI_YELLOW = "#F0E442"
OI_BLACK = "#1A1A1A"
OI_GRAY = "#6B7280"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "svg.fonttype": "none",
})


def save_both(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_fig1_architecture() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0, 4.2)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc, ec, lw=1.4):
        rect = mpatches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.08,rounding_size=0.12",
            facecolor=fc, edgecolor=ec, linewidth=lw,
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h - 0.24, title, ha="center", va="top",
                fontsize=9.2, fontweight="bold", color=OI_BLACK)
        ax.text(x + w / 2, y + h / 2 - 0.14, body, ha="center", va="center",
                fontsize=7.9, color="#2D3748", linespacing=1.28)

    def arrow(x1, y1, x2, y2, color=OI_GRAY):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>,head_width=0.18,head_length=0.28",
                                    color=color, lw=1.4))

    # Stage 1: Multi-tier seeding
    box(0.2, 2.35, 2.15, 1.45, "1. Multi-Tier Seeding",
        "KP-14 Eisenstein spaced\nLASTZ11 / LASTZ12 (T=1,2)\nHSP-gated RY20, W10\n6-frame 9-group AA5",
        "#EDF6FC", OI_BLUE)

    # Stage 2: LIS & LCB
    box(2.75, 2.35, 2.0, 1.45, "2. Chaining & LCBs",
        "O(M log M) LIS chain\nBidirectional LCBs\nInter-LCB w=9 rescue\nMacro Horizon Bubbles",
        "#EDF6FC", OI_BLUE)

    # Stage 3: BOCPD + Banded MDL
    box(5.15, 2.35, 2.35, 1.45, "3. BOCPD + Banded MDL",
        "Beta-Bernoulli BOCPD\nover KP-17 W-channel\n3-state logplus BandGrid\n+ X-drop flank extension",
        "#FFF7E6", OI_ORANGE)

    # Stage 4: Conformal Abstention
    box(7.9, 2.35, 2.35, 1.45, "4. Conformal Risk Gate",
        "Calibrated floor λ̂ = 10 bits\n(α = 0.07, n_cal = 15)\nKeep ΔL ≥ λ̂ blocks\nEmit Horizon Bubbles",
        "#FDEEE8", OI_VERMILION)

    arrow(2.38, 3.07, 2.72, 3.07)
    arrow(4.78, 3.07, 5.12, 3.07)
    arrow(7.53, 3.07, 7.87, 3.07)

    # Bottom row: E-graph, Sweep-line Merge, Interop
    box(1.1, 0.3, 2.6, 1.35, "5. E-Graph Saturation",
        "Congruence-closed e-graph\ninv(inv(X)) ↔ X, seq(X,X) → dup\nParsimony cost extraction",
        "#EAF7F2", OI_GREEN)

    box(4.15, 0.3, 2.6, 1.35, "6. O(N log N) Merge & Guide",
        "Best-score interval netting\nQuery gap-column threading\nTransitive sister-clade rescue",
        "#EAF7F2", OI_GREEN)

    box(7.2, 0.3, 2.6, 1.35, "7. Lossless Interop",
        "Bidirectional GFA ↔ MAF lens\nrGFA (SN:Z, SO:i, SR:i:0)\nPure-Rust HDF5 HAL decoder",
        "#F7EEF4", OI_PURPLE)

    arrow(9.07, 2.32, 5.45, 1.70, OI_VERMILION)
    arrow(3.73, 0.97, 4.12, 0.97, OI_GREEN)
    arrow(6.78, 0.97, 7.17, 0.97, OI_GREEN)

    ax.set_title("Horizon Whole-Genome Alignment Architecture: Seeding, Adaptive MDL, Conformal Abstention, and Graph Merge",
                 loc="left", pad=8, fontsize=10.5)
    save_both(fig, "fig1_architecture")


def make_fig2_benchmarks() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.3), gridspec_kw={"width_ratios": [1.35, 1.0]})

    # Panel A: Recall vs Precision across external suites
    # Format: (label, recall, precision, regime_color, dx, dy)
    points = [
        ("470-way primates (n=6)", 83.7, 98.4, OI_BLUE, -28, -12),
        ("Pf × Pr Chr1 / WG", 85.6, 98.3, OI_BLUE, -32, 6),
        ("100-way chrM (panTro4)", 88.9, 99.2, OI_BLUE, -45, -14),
        ("sars-cov-2.hal (SARS-1)", 76.8, 98.4, OI_BLUE, -55, 6),
        ("cactus447way primates", 74.4, 93.1, OI_SKY, -35, 7),
        ("27_sauropsids birds (n=4)", 71.0, 83.2, OI_GREEN, 6, -4),
        ("92_mammals EPO (n=6)", 67.2, 90.2, OI_GREEN, -42, -13),
        ("100-way chrM (galGal4)", 65.6, 88.6, OI_GREEN, -42, 6),
        ("Pf × Pv Chr1 (--top-k 4)", 50.6, 79.5, OI_ORANGE, 6, -4),
        ("sars-cov-2.hal (MERS/OC43)", 40.7, 67.9, OI_VERMILION, 6, -4),
        ("cactus447way Laurasiatheria", 37.1, 75.9, OI_VERMILION, 6, 5),
        ("hg38 5-way promoters", 22.6, 49.3, OI_PURPLE, 6, 3),
    ]

    for label, rec, prec, col, dx, dy in points:
        ax1.scatter(rec, prec, color=col, s=58, edgecolor="white", linewidth=0.8, zorder=4)
        ax1.annotate(label, (rec, prec), xytext=(dx, dy), textcoords="offset points",
                     fontsize=7.6, color=OI_BLACK)

    ax1.set_xlim(15, 100)
    ax1.set_ylim(40, 104)
    ax1.set_xlabel("Base-Pair Recall vs. Reference Alignment (%)")
    ax1.set_ylabel("Base-Pair Precision (%)")
    ax1.set_title("A. De Novo Pairwise Accuracy Across 7 External Suites", loc="left")
    ax1.grid(True, linestyle=":", alpha=0.35, zorder=0)

    legend_handles = [
        mpatches.Patch(color=OI_BLUE, label="Close (<15% div.)"),
        mpatches.Patch(color=OI_SKY, label="Primates (cactus447way)"),
        mpatches.Patch(color=OI_GREEN, label="Medium (Mammals / Birds / chrM)"),
        mpatches.Patch(color=OI_ORANGE, label="Cross-clade AT/GC (Pf × Pv)"),
        mpatches.Patch(color=OI_VERMILION, label="Distant (Laurasiatheria / CoVs)"),
        mpatches.Patch(color=OI_PURPLE, label="Ultra-short 1 kb promoters"),
    ]
    ax1.legend(handles=legend_handles, loc="lower right", frameon=True, framealpha=0.92,
               fontsize=7.2, edgecolor="#CBD5E0")

    # Panel B: horizon_merge multi-species scaling
    suites = [
        "27_sauropsids\n(23 sp, 50 blk)",
        "92_mammals\n(92 sp, 30 blk)",
        "100-way chrM\n(92 sp, 459 blk)",
        "470-way\n(317 sp, 200 blk)",
        "cactus447way\n(406 sp, 200 blk)",
    ]
    times = [0.62, 1.38, 3.07, 2.53, 3.84]
    y_pos = np.arange(len(suites))
    bars = ax2.barh(y_pos, times, height=0.56, color=OI_BLUE, alpha=0.88, edgecolor=OI_BLACK, linewidth=0.6)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(suites, fontsize=8.0)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 4.7)
    ax2.set_xlabel("Wall-Clock Merge Runtime (s)")
    ax2.set_title("B. Exact O(N log N) Merge Scaling (Up to 406 Species)", loc="left")
    for bar, t in zip(bars, times):
        ax2.text(t + 0.08, bar.get_y() + bar.get_height() / 2, f"{t:.2f} s",
                 va="center", ha="left", fontsize=8.2, fontweight="bold", color=OI_BLACK)

    plt.tight_layout(w_pad=2.4)
    save_both(fig, "fig2_benchmarks")


def make_fig3_abstention_speedup() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.8))

    # Panel A: Epistemic abstention & out-of-domain preservation
    categories = [
        "Cannabis cs10×Cherry_Pie\n(45.5% GC, n=418,663 cols)",
        "Pf×Pv LASTZ MDL-supported\n(n=88,597 pairs)",
        "Pf×Pv LASTZ MDL-rejected\n(n=20,968 pairs)",
    ]
    kept = [99.95, 54.4, 12.3]
    abstained = [0.05, 45.6, 87.7]
    x = np.arange(len(categories))
    w = 0.52

    ax1.bar(x, kept, width=w, color=OI_BLUE, label="Aligned / Retained (%)")
    ax1.bar(x, abstained, width=w, bottom=kept, color=OI_VERMILION, alpha=0.85, label="Abstained / Horizon Bubble (%)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(categories, fontsize=7.8)
    ax1.set_ylim(0, 112)
    ax1.set_ylabel("Share of Reference Columns / Pairs (%)")
    ax1.set_title("A. Calibrated Abstention vs. True-Alignment Retention", loc="left")
    for i, (k, a) in enumerate(zip(kept, abstained)):
        ax1.text(i, k / 2, f"{k:.1f}%" if k < 99 else f"{k:.2f}%", ha="center", va="center",
                 color="white", fontweight="bold", fontsize=8.2)
        if a > 10:
            ax1.text(i, k + a / 2, f"{a:.1f}%", ha="center", va="center",
                     color="white", fontweight="bold", fontsize=8.2)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.5, edgecolor="#CBD5E0")

    # Panel B: Whole-genome runtime speedup on slowest Plasmodium pair (Pm x Pb)
    pipelines = ["Galaxy Growler\n(KegAlign + LASTZ)\n[859.9 min]", "horizon_wga\n(14 × 45 contigs)\n[17.12 min]"]
    mins = [859.9, 17.12]
    colors = [OI_GRAY, OI_GREEN]
    bars = ax2.bar(np.arange(2), mins, width=0.46, color=colors, edgecolor=OI_BLACK, linewidth=0.7)
    ax2.set_xticks(np.arange(2))
    ax2.set_xticklabels(pipelines, fontsize=8.5)
    ax2.set_ylim(0, 980)
    ax2.set_ylabel("Whole-Genome Wall-Clock Runtime (min)")
    ax2.set_title("B. P. malariae × P. brasilianum Speedup (50.2×)", loc="left")
    for bar, m in zip(bars, mins):
        ax2.text(bar.get_x() + bar.get_width() / 2, m + 20, f"{m:.1f} min",
                 ha="center", va="bottom", fontsize=8.8, fontweight="bold", color=OI_BLACK)
    ax2.annotate("50.2× faster\n(89.07M aligned pairs,\n536,041 Horizon Bubbles)",
                 xy=(1, 40), xytext=(0.45, 500),
                 arrowprops=dict(arrowstyle="-|>", color=OI_GREEN, lw=1.4),
                 fontsize=8.2, fontweight="bold", color=OI_GREEN,
                 bbox=dict(boxstyle="round,pad=0.3", fc="#EAF7F2", ec=OI_GREEN, lw=0.9))

    plt.tight_layout(w_pad=2.2)
    save_both(fig, "fig3_abstention_speedup")


if __name__ == "__main__":
    make_fig1_architecture()
    make_fig2_benchmarks()
    make_fig3_abstention_speedup()
    print(f"Generated 3 figures (.png + .svg) in {OUT_DIR}")
