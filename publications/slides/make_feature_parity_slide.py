#!/usr/bin/env python3
"""Generate a slide showing CRISPR/feature matching parity for perturb-seq datasets."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 14,
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "axes.edgecolor": "#cccccc",
    "axes.grid": True,
    "grid.alpha": 0.25,
})

datasets = [
    "A375\n(11 guides, minUMI 10)",
    "UCSF EBs2_2\n(548 guides, minUMI 3)",
    "MSK 30polyKO\n(30 guides, minUMI 2)",
]

crispr_match = [100.0, 98.9, 99.4]
umi_pearson  = [1.000, 0.999, 1.000]

counts_text = [
    "1,083 / 1,083",
    "11,902 / 12,038",
    "23,210 / 23,341",
]

CRISPR_COLOR = "#DC2626"
UMI_COLOR    = "#7C3AED"

x = np.arange(len(datasets))
w = 0.30
offsets = [-(w/2 + 0.02), (w/2 + 0.02)]

fig, ax = plt.subplots(figsize=(16, 9))

bars_c = ax.bar(x + offsets[0], crispr_match, w, label="CRISPR Set-Match (%)",
                color=CRISPR_COLOR, edgecolor="white", linewidth=1.2)
bars_u = ax.bar(x + offsets[1], [v * 100 for v in umi_pearson], w,
                label="UMI Pearson (×100)",
                color=UMI_COLOR, edgecolor="white", linewidth=1.2)

for i, (bar, v, ct) in enumerate(zip(bars_c, crispr_match, counts_text)):
    pct = f"{v:.1f}%" if v < 100 else "100%"
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
            pct, ha="center", va="bottom", fontsize=17, fontweight="bold",
            color="#111827")
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() - 3,
            ct, ha="center", va="top", fontsize=12, color="white", fontstyle="italic")

for bar, v in zip(bars_u, umi_pearson):
    label = f"{v:.3f}" if v < 1.0 else "1.000"
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
            label, ha="center", va="bottom", fontsize=17, fontweight="bold",
            color="#111827")

ax.set_xticks(x)
ax.set_xticklabels(datasets, fontsize=15)
ax.set_ylabel("Percent / Correlation ×100", fontsize=16)
ax.set_ylim(0, 112)
ax.set_title("STAR Suite Feature Calling Parity vs Cell Ranger 9",
             fontsize=24, fontweight="bold", pad=20)
ax.axhline(y=100, color="#9CA3AF", linewidth=0.8, linestyle="--", zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(loc="upper center", fontsize=17, framealpha=0.9, edgecolor="#E5E7EB",
          ncol=2)

fig.tight_layout()

out = "publications/slides/feature_parity_bars.png"
fig.savefig(out, dpi=200, bbox_inches="tight", transparent=True)
print(f"Saved → {out}")

out_pdf = "publications/slides/feature_parity_bars.pdf"
fig.savefig(out_pdf, bbox_inches="tight", transparent=True)
print(f"Saved → {out_pdf}")
