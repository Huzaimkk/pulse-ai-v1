"""
diagram_parts.py — hand-drawn-style parts for the block diagrams (Figs 1, 2)
==========================================================================
Dashed layer frames with a label column, numbered connectors and small line
icons (database, chip, ECG trace, screen, sensors), in the style of an
engineering system diagram. All coordinates are inches (figstyle.canvas).
"""

import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle

import figstyle as S

ICON = S.INK_SOFT
LW = 0.9


def layer(ax, x, y, w, h, label, label_w=0.95):
    """Dashed layer frame with its name in a left-hand column."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.06",
                                fill=False, ec=S.MUTED, lw=0.9, ls=(0, (4, 3)), zorder=0))
    ax.plot([x + label_w, x + label_w], [y + 0.06, y + h - 0.06], color=S.RULE,
            lw=0.8, ls=(0, (2, 2)), zorder=0)
    ax.text(x + label_w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=8, fontweight="bold", color=S.INK, linespacing=1.2)


def conn(ax, x1, y1, x2, y2, label=None, dashed=False, color=S.BLUE, side="right",
         lw=0.9, off=0.05):
    """Numbered connector. Solid = data flow, dashed = validation / control."""
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8,
                                 color=color, lw=lw, ls=(0, (4, 2)) if dashed else "-",
                                 shrinkA=0, shrinkB=0, zorder=1))
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        if abs(x2 - x1) < abs(y2 - y1):          # vertical connector
            ax.text(mx + (off if side == "right" else -off), my, label,
                    ha="left" if side == "right" else "right", va="center",
                    fontsize=S.MIN_PT, color=color, zorder=4,
                    bbox=dict(fc="white", ec="none", pad=0.6))
        else:
            ax.text(mx, my + off, label, ha="center", va="bottom",
                    fontsize=S.MIN_PT, color=color, zorder=4,
                    bbox=dict(fc="white", ec="none", pad=0.6))


# ------------------------------------------------------------------ icons ---
def icon_db(ax, cx, cy, s=0.22, color=ICON):
    w, h, e = s, s * 1.15, s * 0.32
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, fc="white", ec="none", zorder=2))
    ax.plot([cx - w / 2] * 2, [cy - h / 2, cy + h / 2], color=color, lw=LW, zorder=3)
    ax.plot([cx + w / 2] * 2, [cy - h / 2, cy + h / 2], color=color, lw=LW, zorder=3)
    for yy in (cy + h / 2, cy, cy - h / 2):
        th = np.linspace(0, np.pi, 40) if yy != cy + h / 2 else np.linspace(0, 2 * np.pi, 80)
        ax.plot(cx + w / 2 * np.cos(th), yy - e / 2 * np.sin(th), color=color, lw=LW, zorder=3)


def icon_ecg(ax, x, y, w=0.42, h=0.16, color=S.RED):
    t = np.linspace(0, 1, 200)
    sig = (0.12 * np.exp(-((t - 0.18) / 0.035) ** 2)
           - 0.12 * np.exp(-((t - 0.40) / 0.010) ** 2)
           + 1.00 * np.exp(-((t - 0.43) / 0.012) ** 2)
           - 0.25 * np.exp(-((t - 0.46) / 0.012) ** 2)
           + 0.28 * np.exp(-((t - 0.70) / 0.06) ** 2))
    ax.plot(x + w * t, y + h * (sig - sig.min()) / (sig.max() - sig.min()),
            color=color, lw=1.0, zorder=3)


def icon_chip(ax, cx, cy, s=0.26, color=ICON):
    ax.add_patch(Rectangle((cx - s / 2, cy - s / 2), s, s, fc="white", ec=color, lw=LW, zorder=3))
    ax.add_patch(Rectangle((cx - s / 5, cy - s / 5), 2 * s / 5, 2 * s / 5, fc=S.GREY_FILL,
                           ec=color, lw=0.6, zorder=3))
    for k in np.linspace(-s / 3, s / 3, 4):
        for (a, b, c, d) in [(cx + k, cy + s / 2, cx + k, cy + s / 2 + 0.05),
                             (cx + k, cy - s / 2, cx + k, cy - s / 2 - 0.05),
                             (cx + s / 2, cy + k, cx + s / 2 + 0.05, cy + k),
                             (cx - s / 2, cy + k, cx - s / 2 - 0.05, cy + k)]:
            ax.plot([a, c], [b, d], color=color, lw=0.7, zorder=3)


def icon_screen(ax, cx, cy, w=0.34, h=0.22, color=ICON):
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, fc="white", ec=color, lw=LW, zorder=3))
    icon_ecg(ax, cx - w / 2 + 0.04, cy - h / 4, w - 0.08, h / 2, color=S.BLUE)


def icon_electrodes(ax, cx, cy, color=ICON):
    for dx in (-0.10, 0.10):
        ax.add_patch(Circle((cx + dx, cy), 0.07, fc="white", ec=color, lw=LW, zorder=3))
        ax.add_patch(Circle((cx + dx, cy), 0.025, fc=color, ec="none", zorder=3))
    ax.plot([cx - 0.10, cx - 0.10, cx + 0.10, cx + 0.10], [cy + 0.07, cy + 0.15, cy + 0.15, cy + 0.07],
            color=color, lw=0.7, zorder=3)


def icon_optical(ax, cx, cy, color=ICON):
    ax.add_patch(Circle((cx - 0.08, cy), 0.05, fc=S.RED_FILL, ec=S.RED, lw=LW, zorder=3))
    ax.add_patch(Rectangle((cx + 0.03, cy - 0.05), 0.10, 0.10, fc="white", ec=color, lw=LW, zorder=3))
    for a in (60, 90, 120):
        r = np.radians(a)
        ax.plot([cx - 0.08 + 0.07 * np.cos(r), cx - 0.08 + 0.13 * np.cos(r)],
                [cy + 0.07 * np.sin(r), cy + 0.13 * np.sin(r)], color=S.RED, lw=0.7, zorder=3)


def icon_cap(ax, cx, cy, color=ICON):
    th = np.linspace(0, np.pi, 60)
    ax.plot(cx + 0.16 * np.cos(th), cy - 0.05 + 0.16 * np.sin(th), color=color, lw=LW, zorder=3)
    ax.plot([cx - 0.16, cx + 0.16], [cy - 0.05, cy - 0.05], color=color, lw=LW, zorder=3)
    for a in (30, 60, 90, 120, 150):
        r = np.radians(a)
        ax.add_patch(Circle((cx + 0.16 * np.cos(r), cy - 0.05 + 0.16 * np.sin(r)), 0.022,
                            fc=color, ec="none", zorder=3))


def item(ax, x, y, w, h, title, body, icon=None, fill="white", edge=S.MUTED, lw=0.9):
    """A box with an optional icon drawn in a left gutter (text stays clear of it)."""
    S.box(ax, x, y, w, h, "", None, fill=fill, edge=edge, lw=lw)
    g = 0.0
    if icon is not None:
        g = 0.42
        icon(ax, x + g / 2 + 0.02, y + h / 2)
    S.box(ax, x + g, y, w - g, h, title, body, fill="none", edge="none", lw=0)
