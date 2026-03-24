#!/usr/bin/env python3
"""Generate a single-panel bar chart of STAR-suite speedups, colored by assay type."""

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

BULK_COLOR = "#F59E0B"
PERTURB_COLOR = "#2563EB"
FLEX_COLOR = "#10B981"

labels = [
    "Bulk PE\n(no Y-removal)",
    "Bulk PE\n(Y-removal)",
    "A375\n(47M, 2 lib)",
    "UCSF EBs2_2\n(445M, 2 lib)",
    "MSK 30polyKO\n(669M, 3 lib)",
    "Flex SC2300771\n(2M reads)",
]
speedups = [1.7, 2.4, 3.8, 3.2, 6.1, 2.0]
types    = ["bulk", "bulk", "perturb", "perturb", "perturb", "flex"]
color_map = {"bulk": BULK_COLOR, "perturb": PERTURB_COLOR, "flex": FLEX_COLOR}
colors   = [color_map[t] for t in types]

fig, ax = plt.subplots(figsize=(16, 9))

x = np.arange(len(labels))
fill_colors = [c if t != "flex" else "none" for c, t in zip(colors, types)]
edge_colors = [("white" if t != "flex" else FLEX_COLOR) for t in types]
edge_widths = [(1.5 if t != "flex" else 2.5) for t in types]

bars = ax.bar(x, speedups, width=0.82, color=fill_colors, edgecolor=edge_colors,
              linewidth=1.5)
for bar, ew, ec in zip(bars, edge_widths, edge_colors):
    bar.set_linewidth(ew)
    bar.set_edgecolor(ec)

for bar, sp, t in zip(bars, speedups, types):
    label = f"~{sp}x" if t == "flex" else f"{sp}x"
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.12,
            label, ha="center", va="bottom", fontsize=18, fontweight="bold",
            color="#111827")

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=13)
ax.set_ylabel("Speedup (fold)", fontsize=16)
ax.set_ylim(0, 7.5)
ax.set_title("STAR Suite Speedup", fontsize=24, fontweight="bold", pad=20)
ax.axhline(y=1, color="#9CA3AF", linewidth=0.8, linestyle="--", zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=BULK_COLOR, edgecolor="white", label="Bulk RNA-seq  (vs external stepwise)"),
    Patch(facecolor=PERTURB_COLOR, edgecolor="white", label="Perturb-seq  (vs Cell Ranger 9)"),
    Patch(facecolor="none", edgecolor=FLEX_COLOR, linewidth=2, label="Flex  (vs Cell Ranger 7.2, preliminary)"),
]
ax.legend(handles=legend_elements, loc="upper left", fontsize=13, framealpha=0.9,
          edgecolor="#E5E7EB")

fig.tight_layout()

out = "publications/slides/speedup_bars.png"
fig.savefig(out, dpi=200, bbox_inches="tight", transparent=True)
print(f"Saved → {out}")

out_pdf = "publications/slides/speedup_bars.pdf"
fig.savefig(out_pdf, bbox_inches="tight", transparent=True)
print(f"Saved → {out_pdf}")
