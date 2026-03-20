#!/usr/bin/env python3
"""Generate a parity slide showing Jaccard, Gene Pearson, and Cell Pearson."""

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
    "A375\n(47M, 2 lib)",
    "UCSF EBs2_2\n(445M, 2 lib)",
    "MSK 30polyKO\n(669M, 3 lib)",
    "Flex SC2300771\n(full, 4 samples)",
]

jaccard      = [0.976, 0.976, 0.942, 0.98]
gene_pearson = [0.975, 0.995, 0.994, 0.998]
cell_pearson = [0.9995, 1.000, 1.000, 1.000]

JACCARD_COLOR = "#F59E0B"
GENE_COLOR    = "#2563EB"
CELL_COLOR    = "#10B981"

x = np.arange(len(datasets))
n_metrics = 3
w = 0.24
offsets = [-(w + 0.02), 0, (w + 0.02)]

fig, ax = plt.subplots(figsize=(16, 9))

bars_j = ax.bar(x + offsets[0], jaccard,      w, label="Jaccard (cell overlap)",
                color=JACCARD_COLOR, edgecolor="white", linewidth=1.2)
bars_g = ax.bar(x + offsets[1], gene_pearson, w, label="Gene Pearson",
                color=GENE_COLOR, edgecolor="white", linewidth=1.2)
bars_c = ax.bar(x + offsets[2], cell_pearson, w, label="Cell Pearson",
                color=CELL_COLOR, edgecolor="white", linewidth=1.2)

for bars, vals in [(bars_j, jaccard), (bars_g, gene_pearson), (bars_c, cell_pearson)]:
    for bar, v in zip(bars, vals):
        label = f"{v:.3f}" if v < 1.0 else "1.000"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.001,
                label, ha="center", va="bottom", fontsize=16, fontweight="bold",
                color="#374151")

ax.set_xticks(x)
ax.set_xticklabels(datasets, fontsize=14)
ax.set_ylabel("Correlation / Overlap", fontsize=16)
ax.set_ylim(0, 1.12)
ax.set_title("STAR Suite Parity vs Cell Ranger", fontsize=24, fontweight="bold", pad=20)
ax.axhline(y=1.0, color="#9CA3AF", linewidth=0.8, linestyle="--", zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(loc="upper center", fontsize=17, framealpha=0.9, edgecolor="#E5E7EB",
          ncol=3)

fig.tight_layout()

out = "publications/slides/parity_bars.png"
fig.savefig(out, dpi=200, bbox_inches="tight", transparent=True)
print(f"Saved → {out}")

out_pdf = "publications/slides/parity_bars.pdf"
fig.savefig(out_pdf, bbox_inches="tight", transparent=True)
print(f"Saved → {out_pdf}")
