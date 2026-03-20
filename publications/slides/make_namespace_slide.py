#!/usr/bin/env python3
"""Generate a graphic showing TRU-to-NXT barcode namespace translation
with separate GEX and Feature barcode whitelists."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

fig, ax = plt.subplots(figsize=(16, 9))
ax.set_xlim(0, 16)
ax.set_ylim(0, 9)
ax.axis("off")
fig.patch.set_facecolor("none")

FLANK_COLOR = "#DBEAFE"
CENTER_TRU  = "#FDE68A"
CENTER_NXT  = "#A7F3D0"
ARROW_COLOR = "#6B7280"
GEX_ACCENT  = "#2563EB"
FEAT_ACCENT = "#DC2626"
BOX_BG      = "#F9FAFB"

def draw_barcode(ax, x0, y0, bases, center_color, label, label_color, accent):
    """Draw a 16-bp barcode with flanks and highlighted center."""
    bw, bh = 0.55, 0.55
    total_w = 16 * bw
    sx = x0 - total_w / 2

    bg = FancyBboxPatch((sx - 0.15, y0 - 0.15), total_w + 0.3, bh + 0.3,
                        boxstyle="round,pad=0.1", facecolor=BOX_BG,
                        edgecolor=accent, linewidth=2.0, zorder=0)
    ax.add_patch(bg)

    for i, b in enumerate(bases):
        bx = sx + i * bw
        if 7 <= i <= 8:
            fc = center_color
        else:
            fc = FLANK_COLOR
        rect = FancyBboxPatch((bx + 0.03, y0 + 0.03), bw - 0.06, bh - 0.06,
                              boxstyle="round,pad=0.02", facecolor=fc,
                              edgecolor="#9CA3AF", linewidth=0.8, zorder=1)
        ax.add_patch(rect)
        ax.text(bx + bw / 2, y0 + bh / 2, b, ha="center", va="center",
                fontsize=14, fontweight="bold", fontfamily="monospace", zorder=2)

    ax.text(sx - 0.3, y0 + bh / 2, label, ha="right", va="center",
            fontsize=18, fontweight="bold", color=label_color)

    # flank / center annotations
    mid = sx + 8 * bw
    ax.annotate("", xy=(sx, y0 - 0.12), xytext=(sx + 7 * bw, y0 - 0.12),
                arrowprops=dict(arrowstyle="<->", color="#9CA3AF", lw=1))
    ax.text(sx + 3.5 * bw, y0 - 0.28, "Flank 1 (7 bp)", ha="center",
            fontsize=9, color="#6B7280")
    ax.annotate("", xy=(sx + 7 * bw, y0 - 0.12), xytext=(sx + 9 * bw, y0 - 0.12),
                arrowprops=dict(arrowstyle="<->", color="#9CA3AF", lw=1))
    ax.text(mid, y0 - 0.28, "2 bp", ha="center", fontsize=9, color="#6B7280")
    ax.annotate("", xy=(sx + 9 * bw, y0 - 0.12), xytext=(sx + 16 * bw, y0 - 0.12),
                arrowprops=dict(arrowstyle="<->", color="#9CA3AF", lw=1))
    ax.text(sx + 12.5 * bw, y0 - 0.28, "Flank 2 (7 bp)", ha="center",
            fontsize=9, color="#6B7280")

# ── Title ──
ax.text(8, 8.5, "Feature and Expression Libraries Can Use", ha="center", va="center",
        fontsize=28, fontweight="bold", color="#111827")
ax.text(8, 8.0, "Different Barcode Namespaces", ha="center", va="center",
        fontsize=28, fontweight="bold", color="#111827")

# ── GEX library (TRU) — real barcode from 3M-5pgex-jan-2023.txt ──
gex_tru = list("AAACCATTCAAACCAT")   # center (pos 8-9) = T C
draw_barcode(ax, 8, 6.6, gex_tru, CENTER_TRU,
             "GEX (TRU)  ", GEX_ACCENT, GEX_ACCENT)

# ── Feature library (NXT) — same cell, center swapped ──
feat_nxt = list("AAACCATCTAAACCAT")   # center (pos 8-9) = C T
draw_barcode(ax, 8, 4.6, feat_nxt, CENTER_NXT,
             "Feature (NXT)  ", FEAT_ACCENT, FEAT_ACCENT)

# ── Arrow showing same cell, different namespace ──
ax.annotate("", xy=(8, 5.45), xytext=(8, 6.35),
            arrowprops=dict(arrowstyle="<->, head_width=0.25, head_length=0.12",
                            color=ARROW_COLOR, lw=2))
ax.text(8.6, 5.9, "Same cell,\ndifferent namespace", ha="left", va="center",
        fontsize=12, color=ARROW_COLOR, fontstyle="italic")

# ── Key points ──
points_y = 1.2
points_bg = FancyBboxPatch((1.5, points_y - 0.4), 13, 2.6,
                           boxstyle="round,pad=0.2", facecolor="#F3F4F6",
                           edgecolor="#D1D5DB", linewidth=1.5, zorder=0)
ax.add_patch(points_bg)

bullets = [
    ("Feature cell barcode namespace depends on chemistry and assay", FEAT_ACCENT),
    ("GEX cell barcodes are canonical TRU", GEX_ACCENT),
    ("No public master list maps every 10x feature assay to NXT vs TRU", FEAT_ACCENT),
    ("Outputs are translated / canonicalized to TRU by default", GEX_ACCENT),
    ("STAR Suite autodetects and translates if necessary", "#059669"),
]
for i, (bullet, marker_color) in enumerate(bullets):
    ax.text(2.4, points_y + 1.8 - i * 0.48, "\u25CF", ha="center", va="center",
            fontsize=21, color=marker_color)
    ax.text(2.9, points_y + 1.8 - i * 0.48, bullet, ha="left", va="center",
            fontsize=21, color="#1F2937")

fig.tight_layout()

out = "publications/slides/namespace_translation.png"
fig.savefig(out, dpi=200, bbox_inches="tight", transparent=True)
print(f"Saved → {out}")

out_pdf = "publications/slides/namespace_translation.pdf"
fig.savefig(out_pdf, bbox_inches="tight", transparent=True)
print(f"Saved → {out_pdf}")
