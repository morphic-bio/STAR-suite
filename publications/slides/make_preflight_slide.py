#!/usr/bin/env python3
"""Generate a slide about the Perturb-seq FASTQ pre-flight pipeline."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(16, 9))
ax.set_xlim(0, 16)
ax.set_ylim(0, 9)
ax.axis("off")
fig.patch.set_facecolor("none")

GEX_COLOR  = "#2563EB"
FEAT_COLOR = "#DC2626"
JOIN_COLOR = "#7C3AED"
OK_COLOR   = "#059669"
BOX_BG     = "#F9FAFB"
ARROW_CLR  = "#6B7280"

# ── Title ──
ax.text(8, 8.55, "Perturb-seq FASTQ Pre-Flight", ha="center", va="center",
        fontsize=30, fontweight="bold", color="#111827")
ax.text(8, 8.0, "Automated library detection, namespace translation, and lane pairing",
        ha="center", va="center", fontsize=20, color="#6B7280")

# ── Pipeline steps (left to right flow) ──
step_boxes = [
    {"x": 2.2, "label": "Detect\nLibrary Type", "color": GEX_COLOR,
     "sub": "GEX or Feature?"},
    {"x": 6.0, "label": "Translate to\nCommon Space", "color": FEAT_COLOR,
     "sub": "NXT → TRU auto"},
    {"x": 9.8, "label": "Pair Lanes\n& Features", "color": JOIN_COLOR,
     "sub": "Top-500 BC similarity"},
    {"x": 13.6, "label": "Validated\nConfig", "color": OK_COLOR,
     "sub": "Ready to run"},
]

box_w, box_h = 2.8, 1.5
for step in step_boxes:
    bx = step["x"] - box_w / 2
    by = 5.6
    bg = FancyBboxPatch((bx, by), box_w, box_h,
                        boxstyle="round,pad=0.15", facecolor=step["color"],
                        edgecolor="white", linewidth=2, alpha=0.15, zorder=0)
    ax.add_patch(bg)
    border = FancyBboxPatch((bx, by), box_w, box_h,
                            boxstyle="round,pad=0.15", facecolor="none",
                            edgecolor=step["color"], linewidth=2.5, zorder=1)
    ax.add_patch(border)
    ax.text(step["x"], by + box_h / 2 + 0.15, step["label"], ha="center", va="center",
            fontsize=16, fontweight="bold", color=step["color"], zorder=2)
    ax.text(step["x"], by + 0.25, step["sub"], ha="center", va="center",
            fontsize=14, color="#111827", fontstyle="italic", zorder=2)

# ── Arrows between steps ──
for i in range(len(step_boxes) - 1):
    x1 = step_boxes[i]["x"] + box_w / 2 + 0.05
    x2 = step_boxes[i + 1]["x"] - box_w / 2 - 0.05
    ax.annotate("", xy=(x2, 6.35), xytext=(x1, 6.35),
                arrowprops=dict(arrowstyle="->,head_width=0.25,head_length=0.15",
                                color=ARROW_CLR, lw=2.5))

# ── Key points box ──
points_y = 0.6
points_bg = FancyBboxPatch((1.2, points_y - 0.2), 13.6, 3.8,
                           boxstyle="round,pad=0.2", facecolor="#F3F4F6",
                           edgecolor="#D1D5DB", linewidth=1.5, zorder=0)
ax.add_patch(points_bg)

bullets = [
    ("Reads barcode structure to classify each FASTQ as Expression or Feature", GEX_COLOR),
    ("Auto-translates Feature barcodes to GEX namespace (NXT → TRU)", FEAT_COLOR),
    ("Computes top-500 barcode similarity to join lanes and features to GEX", JOIN_COLOR),
    ("Catches and corrects user mixups, mislabels, and omissions", FEAT_COLOR),
    ("Tested on UCSF perturb-seq set — runs in seconds", OK_COLOR),
]
for i, (bullet, color) in enumerate(bullets):
    ax.text(2.0, points_y + 3.15 - i * 0.65, "\u25CF", ha="center", va="center",
            fontsize=22, color=color)
    ax.text(2.5, points_y + 3.15 - i * 0.65, bullet, ha="left", va="center",
            fontsize=21, color="#1F2937")

fig.tight_layout()

out = "publications/slides/preflight_slide.png"
fig.savefig(out, dpi=200, bbox_inches="tight", transparent=True)
print(f"Saved → {out}")

out_pdf = "publications/slides/preflight_slide.pdf"
fig.savefig(out_pdf, bbox_inches="tight", transparent=True)
print(f"Saved → {out_pdf}")
