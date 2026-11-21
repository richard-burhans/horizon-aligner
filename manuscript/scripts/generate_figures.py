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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "xtick.direction": "out",
    "ytick.direction": "out",
})


def despine(ax: plt.Axes) -> None:
    """Remove top and right spines for clean BioVis L-frame axes."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_both(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
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
        "FlatSeedIndex (12 B/entry)\n16-bit prefix buckets + SDUST\nRayon LASTZ12 & 6-frame AA5\nRY20 / W10 / W9 masks",
        "#EDF6FC", OI_BLUE)

    # Stage 2: LIS & LCB
    box(2.75, 2.35, 2.0, 1.45, "2. Chaining & DAG",
        "O(M log M) LIS + BTreeSet\nBidirectional LCBs\nhorizon_wga --graph (POA)\nMacro Horizon Bubbles",
        "#EDF6FC", OI_BLUE)

    # Stage 3: BOCPD + Banded MDL
    box(5.15, 2.35, 2.35, 1.45, "3. BOCPD + Banded MDL",
        "Beta-Bernoulli BOCPD\n3-state logplus BandLayout\nTopological DAG MDL DP\n+ X-drop flank extension",
        "#FFF7E6", OI_ORANGE)

    # Stage 4: Conformal Abstention
    box(7.9, 2.35, 2.35, 1.45, "4. Conformal Risk Gate",
        "Calibrated floor λ̂ = 10 bits\n(α = 0.07, n_cal = 15)\nKeep ΔL ≥ λ̂ blocks\nEmit Horizon Bubbles",
        "#FDEEE8", OI_VERMILION)

    arrow(2.38, 3.07, 2.72, 3.07)
    arrow(4.78, 3.07, 5.12, 3.07)
    arrow(7.53, 3.07, 7.87, 3.07)

    # Bottom row: E-graph, Sweep-line Merge, Pangenome Suite, Interop
    box(0.2, 0.3, 2.25, 1.35, "5. MdlCostModel E-Graph",
        "horizon_resolve / --resolve-dups\ninv(inv(X)) ↔ X, seq(X,X) → dup\nMDL bit-cost extraction",
        "#EAF7F2", OI_GREEN)

    box(2.8, 0.3, 2.35, 1.35, "6. Merge & Progressive",
        "horizon_merge O(N log N) sweep\nhorizon_guide transitive rescue\nAncestral consensus & netting",
        "#EAF7F2", OI_GREEN)

    box(5.5, 0.3, 2.35, 1.35, "7. Streaming GFA Suite",
        "horizon_gfa (stats / extract)\nhorizon_lift (BED/GFF3 liftover)\nhorizon_clean (BOCPD L-cuts)",
        "#F7EEF4", OI_PURPLE)

    box(8.2, 0.3, 2.05, 1.35, "8. Lossless Interop",
        "horizon <subcommand>\nMAF v1 ↔ rGFA v1.1 lens\nPure-Rust HDF5 v0/v1 HAL",
        "#F7EEF4", OI_PURPLE)

    arrow(9.07, 2.32, 1.32, 1.70, OI_VERMILION)
    arrow(2.48, 0.97, 2.77, 0.97, OI_GREEN)
    arrow(5.18, 0.97, 5.47, 0.97, OI_GREEN)
    arrow(7.88, 0.97, 8.17, 0.97, OI_PURPLE)

    ax.set_title("Horizon Suite Architecture: Seeding, Adaptive MDL, Conformal Abstention, E-Graphs, and Streaming Pangenome Tools",
                 loc="left", pad=8, fontsize=10.2)
    save_both(fig, "fig1_architecture")


def make_fig2_benchmarks() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.7), gridspec_kw={"width_ratios": [1.45, 1.0]})

    # Panel A: Recall vs Precision across external suites
    # Format: (label, recall, precision, regime_color, dx, dy, ha)
    points = [
        ("470-way primates (n=6)", 83.7, 98.4, OI_BLUE, -82, -14, "right"),
        ("Pf × Pr Chr1 / WG", 85.6, 98.3, OI_BLUE, -10, 16, "center"),
        ("100-way chrM (panTro4)", 88.9, 99.2, OI_BLUE, 10, -14, "left"),
        ("60_amniotes (panTro)", 87.8, 95.1, OI_BLUE, 10, -24, "left"),
        ("sars-cov-2.hal (SARS-1)", 76.8, 98.4, OI_BLUE, -14, 12, "right"),
        ("cactus447way primates", 74.4, 93.1, OI_SKY, 12, -4, "left"),
        ("92_mammals EPO (n=6)", 67.2, 90.2, OI_GREEN, -14, 6, "right"),
        ("65_fish (orySin)", 60.3, 89.6, OI_GREEN, -14, -2, "right"),
        ("100-way chrM (galGal4)", 65.6, 88.6, OI_GREEN, -14, -14, "right"),
        ("27_sauropsids birds (n=4)", 71.0, 83.2, OI_GREEN, 10, -4, "left"),
        ("Pf × Pv Chr1 (--top-k 4)", 50.6, 79.5, OI_ORANGE, 10, -2, "left"),
        ("sars-cov-2.hal (MERS/OC43)", 40.7, 67.9, OI_VERMILION, 10, -4, "left"),
        ("cactus447way Laurasiatheria", 37.1, 75.9, OI_VERMILION, 10, 4, "left"),
        ("60_amniotes (galGal)", 25.7, 54.2, OI_VERMILION, 10, 8, "left"),
        ("hg38 5-way promoters", 22.6, 49.3, OI_PURPLE, 10, -8, "left"),
    ]

    for label, rec, prec, col, dx, dy, ha in points:
        ax1.scatter(rec, prec, color=col, s=56, edgecolor="white", linewidth=0.8, zorder=4)
        ax1.annotate(
            label,
            (rec, prec),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=7.5,
            color=OI_BLACK,
            arrowprops=dict(arrowstyle="-", color="#9CA3AF", lw=0.65, shrinkA=2, shrinkB=4),
        )

    ax1.set_xlim(10, 108)
    ax1.set_ylim(42, 108)
    ax1.set_xlabel("Base-Pair Recall vs. Reference Alignment (%)")
    ax1.set_ylabel("Base-Pair Precision (%)")
    ax1.set_title("A. De Novo Pairwise Accuracy Across 9 External Suites", loc="left")
    ax1.grid(True, linestyle=":", alpha=0.35, zorder=0)

    legend_handles = [
        mpatches.Patch(color=OI_BLUE, label="Close (<15% div.)"),
        mpatches.Patch(color=OI_SKY, label="Primates (cactus447way)"),
        mpatches.Patch(color=OI_GREEN, label="Medium (Mammals / Birds / chrM)"),
        mpatches.Patch(color=OI_ORANGE, label="Cross-clade AT/GC (Pf × Pv)"),
        mpatches.Patch(color=OI_VERMILION, label="Distant (Laurasiatheria / CoVs)"),
        mpatches.Patch(color=OI_PURPLE, label="Ultra-short 1 kb promoters"),
    ]
    ax1.legend(handles=legend_handles, loc="lower right", frameon=True, framealpha=0.94,
               fontsize=7.3, edgecolor="#CBD5E0")

    # Panel B: horizon_merge multi-species scaling
    suites = [
        "27_sauropsids\n(23 sp, 50 blk)",
        "92_mammals\n(92 sp, 30 blk)",
        "470-way\n(317 sp, 200 blk)",
        "100-way chrM\n(92 sp, 459 blk)",
        "65_fish\n(58 sp, 40 blk)",
        "cactus447way\n(406 sp, 200 blk)",
        "60_amniotes\n(60 sp, 40 blk)",
    ]
    times = [0.62, 1.38, 2.53, 3.07, 3.66, 3.84, 16.61]
    y_pos = np.arange(len(suites))
    bars = ax2.barh(y_pos, times, height=0.56, color=OI_BLUE, alpha=0.88, edgecolor=OI_BLACK, linewidth=0.6)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(suites, fontsize=8.0)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 19.8)
    ax2.set_xlabel("Wall-Clock Merge Runtime (s)")
    ax2.set_title("B. Exact O(N log N) Merge Scaling (Up to 406 Species)", loc="left")
    for bar, t in zip(bars, times):
        ax2.text(t + 0.08, bar.get_y() + bar.get_height() / 2, f"{t:.2f} s",
                 va="center", ha="left", fontsize=8.2, fontweight="bold", color=OI_BLACK)

    plt.tight_layout(w_pad=2.4)
    save_both(fig, "fig2_benchmarks")


def make_fig3_abstention_speedup() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.1))

    # Panel A: Epistemic abstention & out-of-domain preservation
    categories = [
        "Cannabis\ncs10 × Cherry_Pie\n(45.5% GC, n=418,663)",
        "Pf × Pv LASTZ\nMDL-supported\n(n=88,597 pairs)",
        "Pf × Pv LASTZ\nMDL-rejected\n(n=20,968 pairs)",
    ]
    kept = [99.95, 54.4, 12.3]
    abstained = [0.05, 45.6, 87.7]
    x = np.arange(len(categories))
    w = 0.50

    ax1.bar(x, kept, width=w, color=OI_BLUE, label="Aligned / Retained (%)")
    ax1.bar(x, abstained, width=w, bottom=kept, color=OI_VERMILION, alpha=0.85, label="Abstained / Horizon Bubble (%)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(categories, fontsize=7.9)
    ax1.set_ylim(0, 122)
    ax1.set_ylabel("Share of Reference Columns / Pairs (%)")
    ax1.set_title("A. Calibrated Abstention vs. True-Alignment Retention", loc="left")
    for i, (k, a) in enumerate(zip(kept, abstained)):
        ax1.text(i, k / 2, f"{k:.1f}%" if k < 99 else f"{k:.2f}%", ha="center", va="center",
                 color="white", fontweight="bold", fontsize=8.3)
        if a > 10:
            ax1.text(i, k + a / 2, f"{a:.1f}%", ha="center", va="center",
                     color="white", fontweight="bold", fontsize=8.3)
    ax1.legend(loc="upper center", ncol=2, frameon=True, framealpha=0.94,
               fontsize=7.5, edgecolor="#CBD5E0")

    # Panel B: Whole-genome runtime speedup on slowest Plasmodium pair (Pm x Pb)
    pipelines = ["Galaxy Growler\n(KegAlign + LASTZ)\n[859.9 min]", "horizon_wga\n(14 × 45 contigs)\n[17.12 min]"]
    mins = [859.9, 17.12]
    colors = [OI_GRAY, OI_GREEN]
    bars = ax2.bar(np.arange(2), mins, width=0.44, color=colors, edgecolor=OI_BLACK, linewidth=0.7)
    ax2.set_xticks(np.arange(2))
    ax2.set_xticklabels(pipelines, fontsize=8.5)
    ax2.set_ylim(0, 980)
    ax2.set_ylabel("Whole-Genome Wall-Clock Runtime (min)")
    ax2.set_title("B. P. malariae × P. brasilianum Speedup (50.2×)", loc="left")
    for bar, m in zip(bars, mins):
        ax2.text(bar.get_x() + bar.get_width() / 2, m + 18, f"{m:.1f} min",
                 ha="center", va="bottom", fontsize=8.8, fontweight="bold", color=OI_BLACK)
    ax2.annotate("50.2× faster\n(89.07M aligned pairs,\n536,041 Horizon Bubbles)",
                 xy=(0.94, 72), xytext=(0.38, 490),
                 arrowprops=dict(arrowstyle="-|>", color=OI_GREEN, lw=1.4),
                 fontsize=8.2, fontweight="bold", color=OI_GREEN,
                 bbox=dict(boxstyle="round,pad=0.32", fc="#EAF7F2", ec=OI_GREEN, lw=0.9))

    plt.tight_layout(w_pad=2.5)
    save_both(fig, "fig3_abstention_speedup")


def make_fig4_ablation_profile() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4), gridspec_kw={"width_ratios": [1.1, 1.15]})

    # Panel A: Before vs After Search-Space MDL Gating (60_amniotes ~3 Mb)
    species = ["pan_troglodytes\n(Close)", "mus_musculus\n(Distant)", "gallus_gallus\n(Outgroup)"]
    prec_before = [91.20, 17.69, 6.40]
    prec_after = [95.13, 38.53, 54.24]
    time_before = [223.56, 138.59, 93.71]
    time_after = [25.49, 47.32, 26.69]

    x = np.arange(len(species))
    w = 0.34
    b1 = ax1.bar(x - w / 2, prec_before, width=w, color=OI_GRAY, alpha=0.85, edgecolor=OI_BLACK, linewidth=0.6, label="Fixed-threshold baseline")
    b2 = ax1.bar(x + w / 2, prec_after, width=w, color=OI_BLUE, alpha=0.9, edgecolor=OI_BLACK, linewidth=0.6, label="Search-space MDL + gap-only seeding")

    ax1.set_xticks(x)
    ax1.set_xticklabels(species, fontsize=8.3)
    ax1.set_ylim(0, 118)
    ax1.set_ylabel("Base-Pair Precision (%)")
    ax1.set_title("A. Precision & Speedup on ~3 Mb 60_amniotes Slices", loc="left")
    ax1.grid(True, axis="y", linestyle=":", alpha=0.35)
    for idx, (pb, pa, tb, ta) in enumerate(zip(prec_before, prec_after, time_before, time_after)):
        ax1.text(idx - w / 2, pb + 2.0, f"{pb:.1f}%\n({tb:.0f}s)", ha="center", va="bottom", fontsize=7.2, color="#4B5563")
        speedup = tb / ta
        ax1.text(idx + w / 2, pa + 2.0, f"{pa:.1f}%\n({ta:.1f}s, {speedup:.1f}×)", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=OI_BLUE)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.95, fontsize=7.5, edgecolor="#CBD5E0")

    # Panel B: HORIZON_PROFILE stage breakdown (seconds)
    pairs = ["homo_sapiens ×\npan_troglodytes\n(25.5 s total)", "homo_sapiens ×\nmus_musculus\n(47.3 s total)", "homo_sapiens ×\ngallus_gallus\n(26.7 s total)"]
    w12 = np.array([2.57, 2.47, 1.01])
    w10 = np.array([2.46, 1.92, 1.11])
    lz1 = np.array([12.86, 13.77, 6.70])
    multi = np.array([0.00, 22.02, 14.08])
    dp = np.array([1.43, 2.60, 0.69])
    other = np.array([6.17, 4.54, 3.10])

    y = np.arange(len(pairs))
    stages = [
        ("Canonical w=11/12", w12, OI_BLUE),
        ("w=10 fallback", w10, OI_SKY),
        ("LASTZ12 1-trans", lz1, OI_ORANGE),
        ("LASTZ11 + AA5", multi, OI_VERMILION),
        ("BandLayout MDL DP", dp, OI_GREEN),
        ("X-drop + I/O", other, OI_GRAY),
    ]
    left = np.zeros(len(pairs))
    for name, vals, col in stages:
        ax2.barh(y, vals, left=left, height=0.52, color=col, edgecolor=OI_BLACK, linewidth=0.45, label=name)
        left += vals

    ax2.set_yticks(y)
    ax2.set_yticklabels(pairs, fontsize=8.0)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 54)
    ax2.set_xlabel("Wall-Clock Stage Time (s)")
    ax2.set_title("B. HORIZON_PROFILE Stage Breakdown (3 Mb × 3 Mb)", loc="left")
    ax2.grid(True, axis="x", linestyle=":", alpha=0.35)
    ax2.legend(bbox_to_anchor=(0.48, -0.18), loc="upper center", ncol=3, frameon=True, framealpha=0.95, fontsize=7.3, edgecolor="#CBD5E0")

    plt.tight_layout(w_pad=2.4)
    save_both(fig, "fig4_ablation_profile")


def make_fig_c1_bocpd_mdl_trace() -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.6, 4.8), sharex=True, gridspec_kw={"height_ratios": [1.0, 1.15]})
    x = np.arange(0, 600)

    # Synthetic 600 bp isochore + hypervariable bubble trace matching BOCPD dynamics
    pw = np.empty_like(x, dtype=float)
    pw[:220] = 0.81 + 0.018 * np.sin(x[:220] / 14.0)
    # Sharp transition at 220..240 into 45% GC (P(W)=0.55), with BOCPD rapid reset
    pw[220:240] = np.linspace(0.81, 0.55, 20)
    pw[240:] = 0.55 + 0.02 * np.cos(x[240:] / 16.0)

    ax1.plot(x, pw * 100.0, color=OI_BLUE, lw=2.0, label="BOCPD posterior P(W) (% A+T)")
    ax1.axhline(52.5, color=OI_GRAY, ls="--", lw=1.2, label="Fixed HOXD70 assumption (52.5% A+T)")
    ax1.axvspan(310, 410, color=OI_VERMILION, alpha=0.14)
    ax1.axvline(220, color=OI_ORANGE, ls=":", lw=1.3)
    ax1.text(225, 62, "Changepoint reset\n(19% → 45% GC)", fontsize=7.8, color=OI_ORANGE, fontweight="bold")
    ax1.set_ylabel("Local A+T Share (%)")
    ax1.set_ylim(45, 92)
    ax1.set_title("A. Beta-Bernoulli BOCPD Tracking Across an Isochore Boundary", loc="left")
    ax1.legend(loc="upper right", frameon=True, framealpha=0.94, fontsize=7.6, edgecolor="#CBD5E0")
    ax1.grid(True, linestyle=":", alpha=0.35)

    # Bottom track: Per-column posterior P(M_ij | X, Y) + Horizon Bubble abstention
    post = np.ones_like(x, dtype=float) * 0.94
    post[:220] = 0.93 + 0.04 * np.cos(x[:220] / 11.0)
    post[220:310] = 0.90 + 0.04 * np.sin(x[220:310] / 9.0)
    post[305:315] = np.linspace(0.91, 0.07, 10)
    post[315:405] = 0.07 + 0.03 * np.sin(x[315:405] / 5.0)
    post[405:415] = np.linspace(0.08, 0.92, 10)
    post[415:] = 0.92 + 0.04 * np.cos(x[415:] / 13.0)

    ax2.plot(x, post, color=OI_GREEN, lw=2.0, label="MDL column posterior P(M_ij | X, Y)")
    ax2.axhline(0.15, color=OI_VERMILION, ls="--", lw=1.2, label="Posterior floor (0.15)")
    ax2.axvspan(310, 410, color=OI_VERMILION, alpha=0.18, label="Abstained Horizon Bubble (ΔL < λ̂)")
    ax2.text(360, 0.56, "Horizon Bubble\n(ΔL = -8.4 bits < λ̂)", ha="center", va="center",
             fontsize=8.0, fontweight="bold", color=OI_VERMILION,
             bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=OI_VERMILION, lw=0.8))
    ax2.set_xlabel("Reference Coordinate (bp)")
    ax2.set_ylabel("Column Posterior P(M_ij)")
    ax2.set_ylim(0.0, 1.06)
    ax2.set_title("B. MDL Forward-Backward Posterior & Conformal Abstention Gate", loc="left")
    ax2.legend(loc="lower right", frameon=True, framealpha=0.94, fontsize=7.6, edgecolor="#CBD5E0")
    ax2.grid(True, linestyle=":", alpha=0.35)

    plt.tight_layout()
    save_both(fig, "fig_c1_bocpd_mdl_trace")


def make_fig_c2_bandgrid_egraph() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.3), gridspec_kw={"width_ratios": [1.0, 1.2]})

    # Left: BandLayout diagonal band vs flat Vec<f64>
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.invert_yaxis()
    ax1.set_xticks(range(11))
    ax1.set_yticks(range(11))
    ax1.set_xlabel("Query coordinate j")
    ax1.set_ylabel("Target coordinate i")
    ax1.set_title("A. BandLayout Diagonal Band (W = 2, N = M, stride = 5)", loc="left")
    ax1.grid(True, linestyle=":", alpha=0.35)

    for i in range(10):
        for j in range(10):
            if abs(j - i) <= 2:
                col = OI_BLUE if j == i else "#D8ECFA"
                alpha = 0.85 if j == i else 0.9
                ax1.add_patch(mpatches.Rectangle((j, i), 1, 1, facecolor=col, edgecolor="#64748B", lw=0.6, alpha=alpha))
                txt_col = "white" if j == i else OI_BLACK
                ax1.text(j + 0.5, i + 0.5, f"{j - i + 2}", ha="center", va="center", fontsize=7.2, color=txt_col)
            else:
                ax1.add_patch(mpatches.Rectangle((j, i), 1, 1, facecolor="#F3F4F6", edgecolor="#E5E7EB", lw=0.4))

    ax1.text(5.0, 9.55, "Flat index = i · (2W + 1) + (j − i + W)", ha="center", va="center",
             fontsize=8.0, fontweight="bold", color=OI_BLACK,
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=OI_BLUE, lw=0.9))

    # Right: E-Graph Duplication & Inversion Collapse
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 6)
    ax2.axis("off")
    ax2.set_title("B. Congruence-Closed E-Graph Paralog Collapse", loc="left")

    def ebox(bx, by, bw, bh, title, body, fc, ec):
        p = mpatches.FancyBboxPatch((bx, by), bw, bh, boxstyle="round,pad=0.1,rounding_size=0.15",
                           facecolor=fc, edgecolor=ec, linewidth=1.3)
        ax2.add_patch(p)
        ax2.text(bx + bw / 2, by + bh - 0.32, title, ha="center", va="top", fontsize=8.3, fontweight="bold", color=OI_BLACK)
        ax2.text(bx + bw / 2, by + 0.42, body, ha="center", va="center", fontsize=7.5, color="#374151")

    ebox(0.3, 4.1, 4.1, 1.4, "Input Synteny Sequence", "seq(inv(inv(A)), A)\n(raw LCB order)", "#FDEEE8", OI_VERMILION)
    ebox(5.4, 4.1, 4.2, 1.4, "Rule 1: Double Inversion", "inv(inv(A)) ↔ A\n⇒ e-class(A)", "#FFF7E6", OI_ORANGE)
    ebox(0.3, 1.8, 4.1, 1.4, "Congruence Merge", "seq(A, A) ≡ seq(inv(inv(A)), A)\n(same e-class)", "#EDF6FC", OI_BLUE)
    ebox(5.4, 1.8, 4.2, 1.4, "Rule 4: Tandem Collapse", "seq(A, A) → dup(A)\n(# EG:Z: parsimonious tag)", "#EAF7F2", OI_GREEN)

    ax2.annotate("", xy=(5.32, 4.8), xytext=(4.48, 4.8), arrowprops=dict(arrowstyle="-|>", color=OI_BLACK, lw=1.3))
    ax2.annotate("", xy=(2.35, 3.28), xytext=(7.4, 4.04), arrowprops=dict(arrowstyle="-|>", color=OI_BLACK, lw=1.3, connectionstyle="arc3,rad=0.12"))
    ax2.annotate("", xy=(5.32, 2.5), xytext=(4.48, 2.5), arrowprops=dict(arrowstyle="-|>", color=OI_GREEN, lw=1.4))
    ax2.text(5.0, 0.7, "Extracted Minimum-Cost Form: dup(A)   (cost = 1 vs. raw cost = 4)",
             ha="center", va="center", fontsize=8.3, fontweight="bold", color=OI_GREEN,
             bbox=dict(boxstyle="round,pad=0.3", fc="#EAF7F2", ec=OI_GREEN, lw=0.9))

    plt.tight_layout(w_pad=2.2)
    save_both(fig, "fig_c2_bandgrid_egraph")


def make_fig_c3_hdf5_hal_seek() -> None:
    fig, ax = plt.subplots(figsize=(10.8, 3.8))
    ax.set_xlim(0, 11.2)
    ax.set_ylim(0, 4.2)
    ax.axis("off")

    def hbox(x, y, w, h, title, body, fc, ec):
        p = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                           facecolor=fc, edgecolor=ec, linewidth=1.3)
        ax.add_patch(p)
        ax.text(x + w / 2, y + h - 0.16, title, ha="center", va="top", fontsize=8.3, fontweight="bold", color=OI_BLACK)
        ax.text(x + w / 2, y + (h - 0.35) / 2, body, ha="center", va="center", fontsize=7.2, color="#374151")

    def harr(x1, y1, x2, y2, col=OI_BLACK):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3))

    hbox(0.2, 2.35, 2.0, 1.35, "1. Superblock", "8B HDF5 magic\noffset/len = 8B\nroot_group_addr", "#EDF6FC", OI_BLUE)
    hbox(2.6, 2.35, 2.35, 1.35, "2. v1 Object Header", "iter_v1_header_messages_seek\nfollows 0x0010 continuations\nreads SymbolTable (0x0011)", "#EDF6FC", OI_BLUE)
    hbox(5.35, 2.35, 2.45, 1.35, "3. Group B-Tree", "TREE (type 0) → SNOD\n40B entries + HEAP names\n(e.g. Anc00..leaf genomes)", "#EAF7F2", OI_GREEN)
    hbox(8.2, 2.35, 2.75, 1.35, "4. Dataset Metadata", "Dataspace (0x0001): dims\nDatatype (0x0003): stride\nFilterPipeline (0x000B): zlib", "#FFF7E6", OI_ORANGE)

    hbox(1.8, 0.3, 3.6, 1.3, "5a. Contiguous (DataLayout class 1)", "Direct Seek + read_at(data_addr, size)\nused by small arrays (e.g. SEQIDX_ARRAY)", "#F7EEF4", OI_PURPLE)
    hbox(6.0, 0.3, 4.4, 1.3, "5b. Chunked + Deflate (DataLayout class 2)", "Walk TREE (type 1) chunk coords\nbounded ZlibDecoder (≤256 MiB/chunk) → HAL decoders", "#F7EEF4", OI_PURPLE)

    harr(2.24, 3.02, 2.56, 3.02, OI_BLUE)
    harr(4.99, 3.02, 5.31, 3.02, OI_GREEN)
    harr(7.84, 3.02, 8.16, 3.02, OI_ORANGE)
    harr(8.9, 2.31, 3.6, 1.66, OI_PURPLE)
    harr(9.6, 2.31, 8.2, 1.66, OI_PURPLE)

    ax.set_title("Pure-Rust Streaming HDF5 / HAL Seek Traversal (Zero libhdf5 Dependency)", loc="left", pad=6, fontsize=10.5)
    save_both(fig, "fig_c3_hdf5_hal_seek")


def make_fig_g1_hoxd70_vs_mdl_gc() -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.1), gridspec_kw={"width_ratios": [1.15, 1.0]})

    gc = np.linspace(0.12, 0.65, 200)
    at = 1.0 - gc
    q_a = q_t = at / 2.0
    q_c = q_g = gc / 2.0

    # Random identical-base collision probability: q_A^2 + q_T^2 + q_C^2 + q_G^2
    p_collide = (q_a**2 + q_t**2 + q_c**2 + q_g**2) * 100.0
    # Per-base null surprise of A/T in bits: -log2(q_A)
    info_at = -np.log2(q_a)

    ax1.axvspan(15, 25, color="#FDEEE8", alpha=0.65, label="Plasmodium regime (18–22% GC)")
    ax1.axvspan(42, 50, color="#EDF6FC", alpha=0.65, label="Vertebrate / Cannabis (45–48% GC)")
    l1 = ax1.plot(gc * 100, p_collide, color=OI_VERMILION, lw=2.5, label="Random base-collision rate (%)")
    ax1.axhline(25.0, color=OI_GRAY, lw=1.0, ls="--", alpha=0.7)
    p_195 = p_collide[np.argmin(np.abs(gc - 0.195))]
    p_475 = p_collide[np.argmin(np.abs(gc - 0.475))]
    ax1.scatter([19.5, 47.5], [p_195, p_475], color=OI_VERMILION, s=42, zorder=5)
    ax1.annotate("34.3% random identity\n(32.4% A/T, HOXD70 false = 0.872)", xy=(19.5, p_195), xytext=(24.5, 35.8),
                 fontsize=7.7, fontweight="bold", color=OI_VERMILION,
                 arrowprops=dict(arrowstyle="-|>", color=OI_VERMILION, lw=1.1))
    ax1.set_xlabel("Genome GC Content (%)")
    ax1.set_ylabel("Random Identical-Base Collision (%)", color=OI_VERMILION)
    ax1.tick_params(axis="y", labelcolor=OI_VERMILION)
    ax1.set_xlim(12, 65)
    ax1.set_ylim(22, 41)
    ax1.set_title("A. Random Collision Surge vs. Adaptive Null Bits", loc="left")
    ax1.grid(linestyle=":", alpha=0.4)

    ax1b = ax1.twinx()
    ax1b.spines["top"].set_visible(False)
    l2 = ax1b.plot(gc * 100, info_at, color=OI_BLUE, lw=2.5, ls="-.", label="MDL A/T self-information (−log₂ q_A bits)")
    ax1b.set_ylabel("A/T Information Content (bits)", color=OI_BLUE)
    ax1b.tick_params(axis="y", labelcolor=OI_BLUE)
    ax1b.set_ylim(1.0, 2.6)

    lines = l1 + l2
    labels = [ln.get_label() for ln in lines]
    ax1.legend(lines, labels, loc="upper right", fontsize=7.5, framealpha=0.92)

    # Panel B: Per-base information reward at 81% AT (q_A=q_T=0.405, q_G=q_C=0.095)
    bases = ["A–A / T–T\n(81% AT null)", "G–G / C–C\n(81% AT null)"]
    x = np.arange(len(bases))
    w = 0.34
    hoxd_norm = [0.91, 1.00]
    mdl_bits = [-np.log2(0.405), -np.log2(0.095)]

    b1 = ax2.bar(x - w / 2, hoxd_norm, width=w, color=OI_GRAY, edgecolor=OI_BLACK, linewidth=0.7, label="HOXD70 (+91 / +100, normalized)")
    b2 = ax2.bar(x + w / 2, mdl_bits, width=w, color=OI_BLUE, edgecolor=OI_BLACK, linewidth=0.7, label="Horizon MDL info (−log₂ q_b bits)")
    for rect, fmt in [(b1[0], "+91"), (b1[1], "+100"), (b2[0], f"{mdl_bits[0]:.2f} b"), (b2[1], f"{mdl_bits[1]:.2f} b")]:
        ax2.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 0.08, fmt,
                 ha="center", va="bottom", fontsize=8.2, fontweight="bold", color=OI_BLACK)
    ax2.set_xticks(x)
    ax2.set_xticklabels(bases, fontsize=8.5)
    ax2.set_ylabel("Match Weight (HOXD70 norm. / Bits)")
    ax2.set_ylim(0, 4.1)
    ax2.set_title("B. Self-Calibrating Match Weights (81% AT)", loc="left")
    ax2.grid(axis="y", linestyle=":", alpha=0.4)
    ax2.legend(loc="upper left", fontsize=7.8, framealpha=0.92)

    plt.tight_layout(w_pad=2.0)
    save_both(fig, "fig_g1_hoxd70_vs_mdl_gc")


def make_fig_g2_hal_binary_layouts() -> None:
    fig, ax = plt.subplots(figsize=(10.8, 4.6))
    ax.set_xlim(0, 12.0)
    ax.set_ylim(0, 5.4)
    ax.axis("off")

    def row_fields(y, title, fields, total_label):
        ax.text(0.15, y + 0.54, title, fontsize=8.8, fontweight="bold", color=OI_BLACK, va="bottom")
        ax.text(11.85, y + 0.54, total_label, fontsize=8.0, fontweight="bold", color=OI_GRAY, va="bottom", ha="right")
        cur_x = 0.15
        total_w = 11.7
        sum_rel = sum(f[0] for f in fields)
        for rel_w, label, fc, ec in fields:
            w = total_w * (rel_w / sum_rel)
            rect = mpatches.FancyBboxPatch((cur_x, y), w - 0.05, 0.46, boxstyle="round,pad=0.02,rounding_size=0.06",
                                           facecolor=fc, edgecolor=ec, linewidth=1.1)
            ax.add_patch(rect)
            ax.text(cur_x + (w - 0.05) / 2, y + 0.23, label, ha="center", va="center", fontsize=7.4, color=OI_BLACK, fontweight="bold")
            cur_x += w

    row_fields(4.3, "DNA_ARRAY (decode_hal_dna)", [
        (1, "high nibble (4b)\nbase 2k", "#EDF6FC", OI_BLUE),
        (1, "low nibble (4b)\nbase 2k+1", "#D8ECFA", OI_BLUE),
        (2.5, "0..3 = a,c,g,t   |   8..11 = A,C,G,T   |   4 = n   |   else N", "#F8FAFC", OI_GRAY),
    ], "4 bits / base (packed u8)")

    row_fields(3.3, "SEQIDX_ARRAY (decode_hal_seq_index)", [
        (1, "start: u64\n[0..8)", "#EAF7F2", OI_GREEN),
        (1, "top_offset: u64\n[8..16)", "#D6F0E5", OI_GREEN),
        (1, "bottom_offset: u64\n[16..24)", "#EAF7F2", OI_GREEN),
    ], "24 bytes / sequence")

    row_fields(2.3, "TOP_ARRAY (decode_hal_top_array)", [
        (1, "genome_idx: i64\n[0..8)", "#FFF7E6", OI_ORANGE),
        (1, "bottom_idx: i64\n[8..16)", "#FDEFD2", OI_ORANGE),
        (1, "paralogy_idx: i64\n[16..24)", "#FFF7E6", OI_ORANGE),
        (1, "parent_idx: i64\n[24..32)", "#FDEFD2", OI_ORANGE),
        (0.55, "rev: u8\n[32]", "#FDEEE8", OI_VERMILION),
    ], "33 bytes / segment")

    row_fields(1.3, "BOTTOM_ARRAY (decode_hal_bottom_array)", [
        (1, "start: i64\n[0..8)", "#F7EEF4", OI_PURPLE),
        (1, "len: u64\n[8..16)", "#EFE0EA", OI_PURPLE),
        (1, "top_parse: i64\n[16..24)", "#F7EEF4", OI_PURPLE),
        (1.6, "child_0..C-1: (child_idx: i64, rev: u8)\n9 bytes × C children", "#FDEEE8", OI_VERMILION),
    ], "24 + 9C bytes / segment")

    row_fields(0.3, "SEQNAME_ARRAY (decode_hal_seq_names)", [
        (3, "fixed-width UTF-8 / ASCII contig name (trimmed at first 0x00 byte)", "#EDF6FC", OI_BLUE),
        (1, "0x00 padding...", "#F8FAFC", OI_GRAY),
    ], "stride bytes / sequence")

    ax.set_title("HAL Little-Endian Binary Record Layouts Decoded by src/hal.rs", loc="left", pad=6, fontsize=10.5)
    save_both(fig, "fig_g2_hal_binary_layouts")


def make_fig_g3_horizon_bubble_anatomy() -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.6, 5.2), gridspec_kw={"height_ratios": [1.05, 1.0]})

    # Panel A: Two Origins of a Horizon Bubble along a Chromosome Pair
    ax1.set_xlim(0, 12.0)
    ax1.set_ylim(0, 3.2)
    ax1.axis("off")
    ax1.set_title("A. Two Paths to Abstention: Micro (Conformal) vs. Macro (Structural) Horizon Bubbles", loc="left", fontsize=10.2)

    # Target and Query tracks
    ax1.text(0.15, 2.35, "Target X", fontsize=8.5, fontweight="bold", color=OI_BLACK, va="center")
    ax1.text(0.15, 0.85, "Query Y", fontsize=8.5, fontweight="bold", color=OI_BLACK, va="center")

    blocks = [
        (1.3, 2.4, "Aligned LCB 1\n(ΔL ≥ λ_eff)", "#EDF6FC", OI_BLUE, True),
        (3.85, 2.1, "Micro Bubble\nΔL < λ_eff or P(M)<0.15", "#FDEEE8", OI_VERMILION, False),
        (6.1, 2.1, "Aligned LCB 2\n(ΔL ≥ λ_eff)", "#EDF6FC", OI_BLUE, True),
        (8.35, 3.3, "Macro Bubble (> 2.5 kb gap)\nNo w=9 two-hit support (DP skipped)", "#FFF7E6", OI_ORANGE, False),
    ]
    for bx, bw, lbl, fc, ec, aligned in blocks:
        for yy in (2.05, 0.55):
            r = mpatches.FancyBboxPatch((bx, yy), bw, 0.58, boxstyle="round,pad=0.04,rounding_size=0.08",
                                        facecolor=fc, edgecolor=ec, linewidth=1.3, linestyle="-" if aligned else "--")
            ax1.add_patch(r)
        if aligned:
            for dx in np.linspace(0.25, bw - 0.25, 6):
                ax1.plot([bx + dx, bx + dx], [1.18, 2.00], color=OI_BLUE, lw=1.0, alpha=0.65)
            ax1.text(bx + bw / 2, 1.59, lbl, ha="center", va="center", fontsize=7.6, fontweight="bold", color=OI_BLUE,
                     bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=OI_BLUE, lw=0.7))
        else:
            ax1.plot([bx + 0.15, bx + bw - 0.15], [1.18, 2.00], color=ec, lw=1.1, ls=":")
            ax1.plot([bx + 0.15, bx + bw - 0.15], [2.00, 1.18], color=ec, lw=1.1, ls=":")
            ax1.text(bx + bw / 2, 1.59, lbl, ha="center", va="center", fontsize=7.4, fontweight="bold", color=ec,
                     bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=ec, lw=0.8))

    # Panel B: GFA v1.1 Graph vs MAF Projection (Bidirectional Lens)
    ax2.set_xlim(0, 12.0)
    ax2.set_ylim(0, 3.0)
    ax2.axis("off")
    ax2.set_title("B. GFA v1.1 Graph Representation & Bidirectional GFA ↔ MAF Lens Preservation", loc="left", fontsize=10.2)

    nodes = [
        (0.5, 1.65, 2.1, "S  blk_0\n(aligned, SR:i:0)", "#EDF6FC", OI_BLUE),
        (3.3, 1.65, 2.3, "S  bub_1  *\nLN:i:85  HB:i:1", "#FDEEE8", OI_VERMILION),
        (6.3, 1.65, 2.1, "S  blk_2\n(aligned, SR:i:0)", "#EDF6FC", OI_BLUE),
        (9.1, 1.65, 2.4, "S  bub_3  *\nLN:i:3100  HB:i:1", "#FFF7E6", OI_ORANGE),
    ]
    for nx, ny, nw, nlbl, fc, ec in nodes:
        p = mpatches.FancyBboxPatch((nx, ny), nw, 0.78, boxstyle="round,pad=0.05,rounding_size=0.1",
                                    facecolor=fc, edgecolor=ec, linewidth=1.3)
        ax2.add_patch(p)
        ax2.text(nx + nw / 2, ny + 0.39, nlbl, ha="center", va="center", fontsize=7.6, fontweight="bold", color=OI_BLACK)

    for x1, x2 in [(2.66, 3.24), (5.66, 6.24), (8.46, 9.04)]:
        ax2.annotate("", xy=(x2, 2.04), xytext=(x1, 2.04), arrowprops=dict(arrowstyle="-|>", color=OI_BLACK, lw=1.3))
        ax2.text((x1 + x2) / 2, 2.22, "L +/+ 0M", ha="center", va="bottom", fontsize=6.7, color=OI_GRAY)

    maf_box = mpatches.FancyBboxPatch((0.5, 0.25), 11.0, 0.85, boxstyle="round,pad=0.06,rounding_size=0.1",
                                      facecolor="#EAF7F2", edgecolor=OI_GREEN, linewidth=1.2)
    ax2.add_patch(maf_box)
    ax2.text(6.0, 0.68,
             "MAF View A = get(C): emits blk_0 and blk_2 alignment blocks only\n"
             "Lens putback(A', C): preserves bub_1 and bub_3 (HB:i:1) unchanged in constant complement (GETPUT / PUTGET)",
             ha="center", va="center", fontsize=7.8, fontweight="bold", color=OI_GREEN)

    ax2.annotate("", xy=(1.55, 1.15), xytext=(1.55, 1.60), arrowprops=dict(arrowstyle="<->", color=OI_GREEN, lw=1.3))
    ax2.annotate("", xy=(7.35, 1.15), xytext=(7.35, 1.60), arrowprops=dict(arrowstyle="<->", color=OI_GREEN, lw=1.3))

    plt.tight_layout(h_pad=1.6)
    save_both(fig, "fig_g3_horizon_bubble_anatomy")


def make_fig5_pangenome_suite() -> None:
    """Two-panel BioVis figure: chr6 stage speedups & 12 GB HPRC suite scaling."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.8, 4.2), gridspec_kw={"width_ratios": [1.1, 1.25], "wspace": 0.34}
    )

    # Panel A: Full GRCh38.chr6 (170.8 Mbp) x CHM13.chr6 (110.7 Mbp) stage breakdown
    configs = [
        "Baseline\n(1e47fe0)",
        "FlatSeedIndex +\nRayon (53130b8)",
        "DP + seed_end +\nw=10 (8cd5db0)",
    ]
    stages = [
        ("w=12 + w=10", [29.78, 28.46, 11.73], OI_BLUE),
        ("LASTZ12 1-trans (lz1)", [588.67, 47.87, 24.54], OI_ORANGE),
        ("HSP gate + LCBs (seed_end)", [184.65, 91.04, 22.88], OI_SKY),
        ("BandLayout DP + X-drop + I/O", [259.35, 283.98, 208.71], OI_GREEN),
    ]
    x = np.arange(len(configs))
    bottom = np.zeros(len(configs))
    for label, vals, col in stages:
        arr = np.array(vals)
        ax1.bar(x, arr, width=0.48, bottom=bottom, color=col, edgecolor="white", linewidth=0.7, label=label)
        bottom += arr

    ax1.text(0, 1062.45 + 22, "1,062.45 s\n(17m 42s)", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=OI_BLACK)
    ax1.text(1, 451.35 + 22, "451.35 s\n(2.35× faster)", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=OI_BLUE)
    ax1.text(2, 267.86 + 22, "267.86 s\n(3.97× faster)", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=OI_GREEN)
    ax1.annotate(
        "lz1: 588.7s → 47.9s → 24.5s\n(24.0× faster)",
        xy=(2.25, 24.0),
        xytext=(1.70, 770),
        ha="center",
        fontsize=7.2,
        fontweight="bold",
        color=OI_VERMILION,
        arrowprops=dict(arrowstyle="->", color=OI_VERMILION, lw=1.0, connectionstyle="arc3,rad=-0.30"),
    )
    ax1.annotate(
        "seed_end: 184.7s → 22.9s (8.1×)\nDP+X-drop+I/O: 284.0s → 208.7s",
        xy=(2.25, 165.0),
        xytext=(1.70, 620),
        ha="center",
        fontsize=7.1,
        fontweight="bold",
        color=OI_BLACK,
        arrowprops=dict(arrowstyle="->", color=OI_BLACK, lw=0.9, connectionstyle="arc3,rad=-0.28"),
    )
    ax1.set_xticks(x)
    ax1.set_xticklabels(configs, fontsize=7.8)
    ax1.set_xlim(-0.55, 2.55)
    ax1.set_ylim(0, 1250)
    ax1.set_ylabel("Wall-Clock Runtime (s)", fontsize=9.0)
    ax1.set_title("Full GRCh38.chr6 × CHM13.chr6 Breakdown", loc="left", fontsize=9.8)
    ax1.legend(frameon=False, fontsize=7.0, loc="upper right")
    despine(ax1)
    ax1.text(-0.16, 1.07, "A", transform=ax1.transAxes, fontsize=11.5, fontweight="bold", va="top")

    # Panel B: 12 GB HPRC v1.0 Graph & Suite Scaling (Wall Time vs Peak RSS)
    suite_runs = [
        ("horizon_gfa stats\n(12 GB HPRC, 81.4M segs)", 196.93, 96.9, OI_PURPLE, (10, -2), "left"),
        ("horizon_gfa extract\n(GRCh38+CHM13 chr22)", 193.89, 474.2, OI_PURPLE, (10, -14), "left"),
        ("horizon_gfa extract\n(GRCh38+CHM13 chr6)", 317.00, 1526.7, OI_PURPLE, (10, 4), "left"),
        ("horizon_lift --to-chrom chr22\n(500 chr22 windows)", 177.55, 708.8, OI_GREEN, (-12, -16), "right"),
        ("horizon_wga --resolve-dups\n(chr22 50.8 × 35.6 Mbp)", 69.10, 3971.6, OI_BLUE, (10, 12), "left"),
        ("horizon_wga --resolve-dups\n(chr6 8cd5db0 170.8 × 110.7 Mbp)", 267.86, 12533.0, OI_BLUE, (-12, 0), "right"),
        ("horizon_wga --resolve-dups\n(chr6 20 Mbp q-arm)", 12.84, 2409.7, OI_SKY, (12, -6), "left"),
    ]
    for label, wt, rss, col, (dx, dy), ha in suite_runs:
        ax2.scatter(wt, rss, color=col, s=64, edgecolor="white", linewidth=0.9, zorder=4)
        ax2.annotate(
            label,
            (wt, rss),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=7.0,
            color=OI_BLACK,
            arrowprops=dict(arrowstyle="-", color="#9CA3AF", lw=0.6, shrinkA=2, shrinkB=4),
        )

    ax2.set_yscale("log")
    ax2.set_xlim(0, 460)
    ax2.set_ylim(50, 26000)
    ax2.set_xlabel("Wall-Clock Runtime (s)", fontsize=9.0)
    ax2.set_ylabel("Peak RSS Memory (MB, log scale)", fontsize=9.0)
    ax2.set_title("12 GB HPRC v1.0 Suite Runtime & Memory", loc="left", fontsize=9.8)
    ax2.grid(True, which="major", linestyle=":", linewidth=0.5, color="#D1D5DB", alpha=0.7)
    despine(ax2)
    ax2.text(-0.16, 1.07, "B", transform=ax2.transAxes, fontsize=11.5, fontweight="bold", va="top")

    save_both(fig, "fig5_pangenome_suite")


def make_fig6_suite_and_vg() -> None:
    """Two-panel BioVis figure: HPRC 10k-seg suite vs vg & BOCPD LCB gate calibration."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.8, 4.2), gridspec_kw={"width_ratios": [1.2, 1.0], "wspace": 0.34}
    )

    # Panel A: 10,000-segment HPRC chr22 subgraph + multi-species suite (Runtime vs RSS)
    tools = [
        ("horizon_gfa stats\n(counter-only)", 0.08, 22.6, OI_PURPLE, (10, -2), "left"),
        ("horizon_clean", 0.13, 33.4, OI_GREEN, (10, 8), "left"),
        ("horizon_wga --graph\n(20k DAG nodes)", 2.98, 153.9, OI_BLUE, (10, 14), "left"),
        ("vg convert -g -p\n(vg v1.60.0)", 1.8292, 76.6, OI_GRAY, (-10, 6), "right"),
        ("vg stats -z -l -N\n(vg v1.60.0, 53.7× slower)", 4.2984, 77.0, OI_VERMILION, (10, -4), "left"),
        ("horizon_guide\n(5-sp HAL stub)", 0.859, 610.1, OI_SKY, (10, 4), "left"),
        ("horizon_merge\n(92-sp chrM)", 21.78, 26.6, OI_GREEN, (-8, -18), "right"),
        ("horizon_resolve\n(e-graph walks)", 5.571, 23.6, OI_ORANGE, (-10, 12), "right"),
    ]
    for label, wt, rss, col, (dx, dy), ha in tools:
        ax1.scatter(wt, rss, color=col, s=68, edgecolor="white", linewidth=0.9, zorder=4)
        ax1.annotate(
            label,
            (wt, rss),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=6.9,
            color=OI_BLACK,
            arrowprops=dict(arrowstyle="-", color="#9CA3AF", lw=0.6, shrinkA=2, shrinkB=4),
        )
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(0.03, 45)
    ax1.set_ylim(12, 1200)
    ax1.set_xlabel("Wall-Clock Runtime (s, log scale)", fontsize=9.0)
    ax1.set_ylabel("Peak RSS Memory (MB, log scale)", fontsize=9.0)
    ax1.set_title("10k-Seg HPRC Subgraph (vs. vg) & Multi-Species Suite", loc="left", fontsize=9.6)
    ax1.grid(True, which="major", linestyle=":", linewidth=0.5, color="#D1D5DB", alpha=0.7)
    despine(ax1)
    ax1.text(-0.15, 1.07, "A", transform=ax1.transAxes, fontsize=11.5, fontweight="bold", va="top")

    # Panel B: R1-HMM BOCPD LCB gate precision gain & 18%-GC null suppression (seed=999)
    cats = ["18%-GC Null\n(÷100, s=999)", "Horse\n(Precision %)", "Macaque\n(Precision %)", "Marmoset\n(Precision %)"]
    before = [180.77, 65.4, 84.6, 86.1]
    after = [0.0, 67.6, 87.0, 87.7]
    x = np.arange(len(cats))
    w = 0.38
    b1 = ax2.bar(x - w / 2, before, width=w, color=OI_GRAY, edgecolor="white", label="Before BOCPD LCB Gate")
    b2 = ax2.bar(x + w / 2, after, width=w, color=OI_BLUE, edgecolor="white", label="After BOCPD LCB Gate (HEAD)")

    for bar, val in zip(b1, before):
        txt = "18,077" if val > 100 else f"{val:.1f}%"
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 2.5, txt, ha="center", va="bottom", fontsize=6.8, color=OI_BLACK)
    for bar, val in zip(b2, after):
        txt = "0" if val == 0 else f"{val:.1f}%"
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 2.5, txt, ha="center", va="bottom", fontsize=6.8, fontweight="bold", color=OI_BLUE)

    ax2.set_xticks(x)
    ax2.set_xticklabels(cats, fontsize=7.5)
    ax2.set_ylim(0, 205)
    ax2.set_ylabel("Precision (%) / Null Pairs (×100)", fontsize=8.8)
    ax2.set_title("BOCPD Composition-Aware LCB Gate (R1-HMM)", loc="left", fontsize=9.6)
    ax2.legend(frameon=False, fontsize=7.2, loc="upper right")
    despine(ax2)
    ax2.text(-0.15, 1.07, "B", transform=ax2.transAxes, fontsize=11.5, fontweight="bold", va="top")

    save_both(fig, "fig6_suite_and_vg")


def make_fig7_chr6_hprc_comparison() -> None:
    """Three-panel BioVis figure: GRCh38.chr6 x CHM13.chr6 Horizon vs minimap2 vs HPRC graph."""
    import json

    data_dir = Path(__file__).resolve().parent.parent / "data"
    hor_npz = data_dir / "chr6_hprc_comparison.npz"
    hor_json = data_dir / "chr6_hprc_comparison.json"
    mm_npz = data_dir / "chr6_minimap2_data.npz"
    mm_json = data_dir / "chr6_minimap2_summary.json"
    if not (hor_npz.exists() and hor_json.exists() and mm_npz.exists() and mm_json.exists()):
        print("Skipping fig7_chr6_hprc_comparison (manuscript/data/chr6_*.npz not found)")
        return

    hdata = np.load(hor_npz)
    mdata = np.load(mm_npz)
    with open(hor_json) as f:
        hsum = json.load(f)
    with open(mm_json) as f:
        msum = json.load(f)

    x_hprc = hdata["x_hprc_bp"] / 1e6
    x_hor_exact = hdata["x_exact_bp"] / 1e6
    x_mm_exact = mdata["x_mm_exact_bp"] / 1e6
    truth_segs = hdata["truth_segs"]
    hor_blocks = hdata["hor_blocks"]

    width_2col = 170.0 * 0.0393701
    fig = plt.figure(figsize=(width_2col, 4.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.15, 1.0], hspace=0.46, wspace=0.34)

    # Panel A: Whole-chromosome synteny dot-plot
    ax_a = fig.add_subplot(gs[0, 0])
    tx_mid = (truth_segs[:, 0] + truth_segs[:, 1]) / 2.0 / 1e6
    ty_mid = ((truth_segs[:, 2] + truth_segs[:, 3]) / 2.0 + hsum["y_offset_bp"]) / 1e6
    ax_a.scatter(tx_mid, ty_mid, s=1.2, c=OI_SKY, alpha=0.35, edgecolors="none", rasterized=True, label="HPRC 1-to-1 segs (≥500 bp)")

    long_mask = (hor_blocks[:, 1] - hor_blocks[:, 0]) >= 2000
    hb_plot = hor_blocks[long_mask] if np.any(long_mask) else hor_blocks
    hx_mid = (hb_plot[:, 0] + hb_plot[:, 1]) / 2.0 / 1e6
    hy_mid = ((hb_plot[:, 2] + hb_plot[:, 3]) / 2.0 + hsum["y_offset_bp"]) / 1e6
    ax_a.scatter(hx_mid, hy_mid, s=1.0, c=OI_BLUE, alpha=0.55, edgecolors="none", rasterized=True, label="Horizon blocks (≥2 kb)")

    ax_a.axvline(61.0, color=OI_ORANGE, linestyle="--", linewidth=0.9, alpha=0.9)
    ax_a.text(63.5, 70.0, "q-arm\nstart", fontsize=6.5, color=OI_VERMILION, va="bottom",
              bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
    ax_a.set_xlim(0, 171)
    ax_a.set_ylim(61, 173)
    ax_a.set_xlabel("GRCh38.chr6 position (Mbp)", fontsize=8.5)
    ax_a.set_ylabel("CHM13.chr6 position (Mbp)", fontsize=8.5)
    ax_a.tick_params(labelsize=7.5)
    ax_a.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=6.4, markerscale=3.0, handletextpad=0.3)
    despine(ax_a)
    ax_a.text(-0.18, 1.08, "A", transform=ax_a.transAxes, fontsize=11, fontweight="bold", va="top")

    # Panel B: Horizon (many-to-many) vs minimap2 -x asm5 (1-to-1 primary) vs HPRC 1-to-1 truth
    ax_b = fig.add_subplot(gs[0, 1])
    metrics = ["Base cov.\nrecall", "Exact (0 bp)\nrecall", "≤ 1 kb\nrecall", "Exact (0 bp)\npair prec."]
    hor_vals = [
        100.0 * hsum["x_base_coverage_recall_of_hprc"],
        100.0 * hsum["exact_coord_recall_on_hprc_1to1"],
        100.0 * hsum["within_1kb_recall_on_hprc_1to1"],
        100.0 * hsum["exact_coord_precision_on_shared"],
    ]
    mm_vals = [
        100.0 * (msum["pairs_on_hprc_1to1_x"] / hsum["hprc_1to1_x_bp"]),
        100.0 * msum["exact_coord_recall_on_hprc_1to1"],
        100.0 * msum["within_1kb_recall_on_hprc_1to1"],
        100.0 * msum["exact_coord_precision_on_shared"],
    ]
    y_pos = np.arange(len(metrics))
    h = 0.34
    b1 = ax_b.barh(y_pos - h / 2, hor_vals, height=h, color=OI_BLUE, alpha=0.92, label="Horizon (many-to-many)")
    b2 = ax_b.barh(y_pos + h / 2, mm_vals, height=h, color=OI_ORANGE, alpha=0.92, label="minimap2 -x asm5 (1-to-1)")
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(metrics, fontsize=7.3)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 118)
    ax_b.set_xlabel("Agreement with HPRC 1-to-1 segs (%)", fontsize=8.5)
    ax_b.tick_params(labelsize=7.5)
    for bar, v in zip(b1, hor_vals):
        ax_b.text(v + 1.2, bar.get_y() + bar.get_height() / 2, f"{v:.1f}%", va="center", ha="left", fontsize=6.4, color=OI_BLUE, fontweight="bold")
    for bar, v in zip(b2, mm_vals):
        ax_b.text(v + 1.2, bar.get_y() + bar.get_height() / 2, f"{v:.1f}%", va="center", ha="left", fontsize=6.4, color=OI_VERMILION, fontweight="bold")
    ax_b.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.01), ncol=1, frameon=False, fontsize=6.4, handlelength=1.4)
    despine(ax_b)
    ax_b.text(-0.22, 1.14, "B", transform=ax_b.transAxes, fontsize=11, fontweight="bold", va="top")

    # Panel C: 1-Mbp binned exact coordinate recall across GRCh38.chr6
    ax_c = fig.add_subplot(gs[1, :])
    bins_mb = np.arange(len(x_hprc)) + 0.5
    ax_c.fill_between(bins_mb, 0, x_hprc, color="#D1D5DB", alpha=0.55, label="HPRC 1-to-1 shared bp")
    ax_c.plot(bins_mb, x_mm_exact, color=OI_ORANGE, linewidth=1.1, linestyle="--", label="minimap2 -x asm5 exact (Δy = 0 bp)")
    ax_c.plot(bins_mb, x_hor_exact, color=OI_BLUE, linewidth=1.15, linestyle="-", label="Horizon exact (Δy = 0 bp)")
    ax_c.axvspan(58.0, 62.0, color=OI_ORANGE, alpha=0.14, label="Centromere / q-arm")
    ax_c.set_xlim(0, 171)
    ax_c.set_ylim(0, 1.05)
    ax_c.set_xlabel("GRCh38.chr6 position (1 Mbp bins, n = 171)", fontsize=8.5)
    ax_c.set_ylabel("Exact match fraction / Mbp", fontsize=8.5)
    ax_c.tick_params(labelsize=7.5)
    ax_c.legend(loc="upper left", bbox_to_anchor=(0.01, 0.96), ncol=1, frameon=False, fontsize=6.4, handlelength=1.5)
    despine(ax_c)
    ax_c.text(-0.08, 1.15, "C", transform=ax_c.transAxes, fontsize=11, fontweight="bold", va="top")

    save_both(fig, "fig7_chr6_hprc_comparison")


def make_fig8_multispecies_scaling() -> None:
    """Two-panel BioVis figure: horizon_guide multi-species scaling & 5-mammal chrM HAL accuracy."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.8, 4.2), gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.34}
    )

    # Panel A: Multi-species horizon_guide Star vs --progressive --hal runtime & RSS
    regimes = [
        "5-sp chrM\n(16.5 kb/sp)",
        "77-sp cactus447way\n(chr22_random, ≥25 kb)",
        "444-sp cactus447way\n(chr22_random, ≥1 kb)",
        "14-sp cactus241way\n(chr22 1MB, 13.2 Mbp)",
    ]
    star_s = [None, 1.7182, 12.8684, 48.0048]
    prog_s = [0.2653, 3.0102, None, 84.8240]
    star_rss = [None, "778 MB", "1.12 GB", "1.80 GB"]
    prog_rss = ["2,109 segs", "803 MB", None, "3.67 GB"]

    y = np.arange(len(regimes))
    h = 0.35
    for i in range(len(regimes)):
        if star_s[i] is not None:
            b = ax1.barh(
                y[i] - h / 2,
                star_s[i],
                height=h,
                color=OI_BLUE,
                edgecolor="white",
                label="Star mode" if i == 1 else None,
            )
            ax1.text(
                star_s[i] * 1.15,
                y[i] - h / 2,
                f"{star_s[i]:.2f} s ({star_rss[i]})",
                va="center",
                ha="left",
                fontsize=6.8,
                color=OI_BLUE,
                fontweight="bold",
            )
        if prog_s[i] is not None:
            pos_y = y[i] if star_s[i] is None else y[i] + h / 2
            b = ax1.barh(
                pos_y,
                prog_s[i],
                height=h,
                color=OI_GREEN,
                edgecolor="white",
                label="--progressive --hal" if i == 0 else None,
            )
            ax1.text(
                prog_s[i] * 1.15,
                pos_y,
                f"{prog_s[i]:.2f} s ({prog_rss[i]})",
                va="center",
                ha="left",
                fontsize=6.8,
                color=OI_GREEN,
                fontweight="bold",
            )

    ax1.set_xscale("log")
    ax1.set_xlim(0.1, 1400)
    ax1.set_yticks(y)
    ax1.set_yticklabels(regimes, fontsize=7.6)
    ax1.invert_yaxis()
    ax1.set_xlabel("Wall-Clock Runtime (s, log scale)", fontsize=9.0)
    ax1.set_title("horizon_guide Scaling (Star vs. Progressive HAL)", loc="left", fontsize=9.2)
    ax1.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax1.grid(True, axis="x", which="major", linestyle=":", linewidth=0.5, color="#D1D5DB", alpha=0.7)
    despine(ax1)
    ax1.text(-0.24, 1.07, "A", transform=ax1.transAxes, fontsize=11.5, fontweight="bold", va="top")

    # Panel B: 5-mammal UCSC chrM.maf hg38-projected progressive HAL accuracy
    taxa = ["panTro4\n(Chimp)", "rheMac3\n(Macaque)", "mm10\n(Mouse)", "canFam3\n(Dog)"]
    recalls = [95.5, 89.2, 71.9, 74.1]
    precs = [99.8, 97.1, 97.5, 98.2]
    x = np.arange(len(taxa))
    w = 0.36
    b1 = ax2.bar(x - w / 2, recalls, width=w, color=OI_SKY, edgecolor="white", label="Recall (%)")
    b2 = ax2.bar(x + w / 2, precs, width=w, color=OI_BLUE, edgecolor="white", label="Precision (%)")

    for bar, val in zip(b1, recalls):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            val + 1.5,
            f"{val:.1f}%",
            ha="center",
            va="bottom",
            fontsize=6.9,
            color=OI_BLACK,
        )
    for bar, val in zip(b2, precs):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            val + 1.5,
            f"{val:.1f}%",
            ha="center",
            va="bottom",
            fontsize=6.9,
            fontweight="bold",
            color=OI_BLUE,
        )

    ax2.set_xticks(x)
    ax2.set_xticklabels(taxa, fontsize=7.8)
    ax2.set_ylim(0, 115)
    ax2.set_ylabel("Pairwise Base-Pair Accuracy (%)", fontsize=9.0)
    ax2.set_title("5-Mammal chrM Progressive HAL (0.265 s, 2,109 segs)", loc="left", fontsize=9.5)
    ax2.legend(frameon=False, fontsize=7.5, loc="lower left")
    despine(ax2)
    ax2.text(-0.15, 1.07, "B", transform=ax2.transAxes, fontsize=11.5, fontweight="bold", va="top")

    save_both(fig, "fig8_multispecies_scaling")


if __name__ == "__main__":
    make_fig1_architecture()
    make_fig2_benchmarks()
    make_fig3_abstention_speedup()
    make_fig4_ablation_profile()
    make_fig5_pangenome_suite()
    make_fig6_suite_and_vg()
    make_fig7_chr6_hprc_comparison()
    make_fig8_multispecies_scaling()
    make_fig_c1_bocpd_mdl_trace()
    make_fig_c2_bandgrid_egraph()
    make_fig_c3_hdf5_hal_seek()
    make_fig_g1_hoxd70_vs_mdl_gc()
    make_fig_g2_hal_binary_layouts()
    make_fig_g3_horizon_bubble_anatomy()
    print(f"Generated 14 figures (.png + .svg + .pdf) in {OUT_DIR}")
