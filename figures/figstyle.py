"""
figstyle.py — one style for every figure in the paper
=====================================================
Import this in every figure script. Nothing else sets a font, a size or a colour.

Why a shared module: IEEE reviewers see the figures as a set. Inconsistent fonts
and box sizes across figures read as carelessness before a single caption is
read. One style file makes consistency automatic rather than a thing to remember.

Design constraints baked in:
  - IEEE two-column: 3.5 in single-column, 7.16 in double-column
  - Minimum MIN_PT (7.5 pt) text at final print size, enforced by check_layout()
  - Vector PDF output, never raster
  - Greyscale-safe: every fill is paired with a distinct EDGE and a TEXT label,
    so the figure survives being printed in black and white

Diagram canvases are laid out in INCHES (canvas() sets the axes limits to the
figure size), and text inside boxes is positioned in POINTS. Both are physical
units, so a font-size change cannot silently push text out of its box: if it
does, check_layout() reports it and save() refuses to write the figure.
"""

from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ---------------------------------------------------------------- geometry ---
COL_1 = 3.5       # inches, IEEE single column
COL_2 = 7.16      # inches, IEEE double column (full width)
MIN_PT = 7.5      # smallest text allowed at final print size

# ------------------------------------------------------------------ colour ---
# From a validated categorical palette. Each is paired with a darker edge so the
# boxes remain distinguishable in greyscale.
INK        = "#0b0b0b"
INK_SOFT   = "#52514e"
MUTED      = "#898781"
RULE       = "#c3c2b7"
SURFACE    = "#ffffff"

BLUE       = "#2a78d6"; BLUE_FILL   = "#e4eefb"
ORANGE     = "#eb6834"; ORANGE_FILL = "#fbe8e0"
AQUA       = "#1baf7a"; AQUA_FILL   = "#dcf2e9"
VIOLET     = "#4a3aa7"; VIOLET_FILL = "#e5e2f3"
RED        = "#d03b3b"; RED_FILL    = "#f9e3e3"
GREY       = "#898781"; GREY_FILL   = "#f0efec"

TITLE_PT = 8.0
BODY_PT = 7.5
LINESPACING = 1.3


def apply():
    """Global rcParams. Call once at the top of every figure script."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        # Arial first: the IEEE figure convention, and narrower than DejaVu Sans,
        # which matters once nothing may be set below MIN_PT.
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "xtick.labelsize": MIN_PT,
        "ytick.labelsize": MIN_PT,
        "legend.fontsize": MIN_PT,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "axes.edgecolor": RULE,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK_SOFT,
        "ytick.color": INK_SOFT,
        "axes.grid": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,      # embed TrueType, required by IEEE
        "ps.fonttype": 42,
    })


# ------------------------------------------------------------------- parts ---
def _text_in(ax, patch, x, y, dy_pt, s, **kw):
    """Text at (x, y) data coords shifted dy_pt points vertically, owned by patch."""
    t = ax.annotate(s, xy=(x, y), xytext=(0, dy_pt), textcoords="offset points",
                    zorder=3, **kw)
    t._container = patch
    return t


def box(ax, x, y, w, h, title, body=None, fill=GREY_FILL, edge=MUTED,
        lw=1.0, title_size=TITLE_PT, body_size=BODY_PT, radius=0.05,
        valign="center", pad_pt=5, gap_pt=2.5, linespacing=LINESPACING):
    """
    A rounded box (x, y, w, h in inches) with a bold title and optional body.

    valign="center" centres the title+body block vertically; "top" anchors it
    pad_pt below the top edge, for tall boxes holding lists.
    """
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=fill, edgecolor=edge, linewidth=lw, zorder=2,
    )
    ax.add_patch(patch)
    cx = x + w / 2

    n_body = body.count("\n") + 1 if body else 0
    body_h = n_body * body_size * linespacing if body else 0.0
    block = title_size + (gap_pt + body_h if body else 0.0)

    if valign == "top":
        anchor_y, top_off = y + h, -pad_pt
    else:
        anchor_y, top_off = y + h / 2, block / 2

    _text_in(ax, patch, cx, anchor_y, top_off, title, ha="center", va="top",
             fontsize=title_size, fontweight="bold", color=INK)
    if body:
        _text_in(ax, patch, cx, anchor_y, top_off - title_size - gap_pt, body,
                 ha="center", va="top", fontsize=body_size, color=INK_SOFT,
                 linespacing=linespacing)
    return patch


def arrow(ax, x1, y1, x2, y2, label=None, color=MUTED, lw=1.0,
          label_size=MIN_PT, label_offset=0.05, style="-|>"):
    """A straight arrow, optionally labelled at its midpoint."""
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle=style, mutation_scale=9,
        color=color, linewidth=lw, zorder=1,
        shrinkA=0, shrinkB=0,
    ))
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        vertical = abs(y2 - y1) > abs(x2 - x1)
        ax.text(mx + (label_offset if vertical else 0),
                my + (0 if vertical else label_offset),
                label, ha="left" if vertical else "center",
                va="center" if vertical else "bottom",
                fontsize=label_size, color=INK_SOFT, zorder=3)


def canvas(width, height):
    """A blank figure whose data coordinates are inches — for diagrams."""
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")
    return fig, ax


# -------------------------------------------------------------- validation ---
def check_layout(fig, min_pt=MIN_PT, tol_px=0.5):
    """
    List every layout problem in a figure; an empty list means it is clean.

    Checks: text below min_pt; text outside the figure (it would widen a
    figure that must be exactly COL_1 or COL_2 wide); text escaping the box it
    belongs to; two pieces of text overlapping.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fig_bb = fig.bbox
    items = []
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_text().strip():
                items.append((t, t.get_window_extent(r)))

    problems = []
    for t, bb in items:
        label = repr(t.get_text().splitlines()[0][:40])
        if t.get_fontsize() < min_pt:
            problems.append(f"{label}: {t.get_fontsize():g} pt < {min_pt:g} pt")
        if (bb.x0 < fig_bb.x0 - tol_px or bb.x1 > fig_bb.x1 + tol_px
                or bb.y0 < fig_bb.y0 - tol_px or bb.y1 > fig_bb.y1 + tol_px):
            problems.append(f"{label}: extends outside the figure")
        box_patch = getattr(t, "_container", None)
        if box_patch is not None:
            cb = box_patch.get_window_extent(r)
            if (bb.x0 < cb.x0 + tol_px or bb.x1 > cb.x1 - tol_px
                    or bb.y0 < cb.y0 + tol_px or bb.y1 > cb.y1 - tol_px):
                problems.append(f"{label}: escapes its box")
    for (t1, b1), (t2, b2) in combinations(items, 2):
        if b1.overlaps(b2):
            problems.append(f"{t1.get_text()[:25]!r} overlaps {t2.get_text()[:25]!r}")
    return problems


def save(fig, name, outdir="fig"):
    """Vector PDF plus a PNG for quick viewing. PDF is what goes in the paper.

    Refuses to save a figure that fails check_layout()."""
    problems = check_layout(fig)
    if problems:
        plt.close(fig)
        raise ValueError(f"{name}: layout problems\n  " + "\n  ".join(problems))
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    pdf, png = out / f"{name}.pdf", out / f"{name}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=200)
    plt.close(fig)
    print(f"  {pdf}  +  {png}")
    return pdf
