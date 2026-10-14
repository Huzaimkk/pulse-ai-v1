"""
make_results_figures.py — results figures of the PULSE AI paper, from saved outputs.

Every number drawn here is read from a CSV written by a re-runnable script in the
repo (out/). Nothing is typed in by hand.

    python make_results_figures.py --out-dir ../out

One colour system across every figure (figstyle colours only):
    RED   = atrial axis (P-wave markers)          -> LAE criteria, pre-stroke axis
    BLUE  = ventricular axis (QRS / QT / ST)      -> arrhythmia axis
    GREY  = reference: controls, rate/demographic markers, non-highlighted models
    INK   = the headline model (all markers)
Within a group comparison, cases take the axis colour and controls are always grey.

Fig 4  cohort flow                         gate1.csv, qc_flow_designA.csv
Fig 5  Bland-Altman, P duration vs LUDB    gate2_beats.csv
Fig 6  case vs control distributions       qc_designA.csv, cohort.csv, group_comparison_designA.csv
Fig 7  effect sizes, all markers           group_comparison_designA.csv
Fig 8  ROC, highlighted best model         models_oof_predictions.csv, models_metrics.csv,
                                           models_permutation_test.csv
Fig 9  permutation importance              models_importance.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
from sklearn.metrics import roc_curve

import figstyle as S

S.apply()

ATRIAL = {"p_duration_ms", "ptfv1_uv_ms", "p_axis_deg", "iab_class", "mvp_score"}
REF = {"mean_hr_bpm", "rmssd_ms", "n_beats", "age", "sex"}


def axis_colour(f):
    if f in ATRIAL:
        return S.RED, S.RED_FILL
    if f in REF:
        return S.GREY, S.GREY_FILL
    return S.BLUE, S.BLUE_FILL


MODEL_LABEL = {"logistic_l2": "Logistic regression", "decision_tree": "Decision tree",
               "random_forest": "Random forest", "hist_gradient_boosting": "Gradient boosting"}
FEATURE_LABEL = {
    "p_duration_ms": "P duration", "ptfv1_uv_ms": "PTFV1", "p_axis_deg": "P axis",
    "mean_hr_bpm": "Mean HR", "rmssd_ms": "RMSSD", "n_beats": "Beats / 10 s",
    "qrs_duration_ms": "QRS duration", "qt_ms": "QT", "qtc_bazett_ms": "QTc (Bazett)",
    "qtc_fridericia_ms": "QTc (Fridericia)", "tpe_ms": "Tpeak–Tend",
    "tpe_qt_ratio": "Tpe/QT", "age": "Age", "sex": "Sex", "mvp_score": "MVP score",
    "iab_class": "IAB class", "fqrs_present": "fQRS",
}


def flabel(f):
    if f.startswith("st_j60_"):
        lead = f[len("st_j60_"):-len("_mv")]
        lead = {"avr": "aVR", "avl": "aVL", "avf": "aVF"}.get(lead, lead.upper())
        return f"ST J+60, {lead}"
    return FEATURE_LABEL.get(f, f)


def clean(ax, left=True):
    for s in ("top", "right") + (() if left else ("left",)):
        ax.spines[s].set_visible(False)


# ===========================================================================
# FIG 4 — cohort flow (canvas units are inches)
# ===========================================================================
def fig04_flow(out):
    g1 = pd.read_csv(out / "gate1.csv").set_index("design").loc["A"]
    fl = pd.read_csv(out / "qc_flow_designA.csv")
    get = lambda grp, rule, col="any": int(fl[(fl.group == grp) & (fl.rule == rule)][col].iloc[0])

    fig, ax = S.canvas(S.COL_1, 3.9)
    x, w = 0.05, 1.55
    cx = x + w / 2
    boxes = [
        (3.38, 0.46, "PTB-XL, 500 Hz", "21,799 records", S.GREY_FILL, S.MUTED),
        (2.58, 0.58, "Eligible",
         f"LAO/LAE {int(g1.cases_eligible):,}\ncontrols {int(g1.controls_eligible):,}",
         S.GREY_FILL, S.MUTED),
        (1.78, 0.58, "Age/sex matched 1:2",
         f"cases {int(g1.cases_matched)}\ncontrols {int(g1.controls_matched)}",
         S.GREY_FILL, S.MUTED),
        (0.05, 0.58, "Analysed",
         f"cases {get('case', 'final_n')}\ncontrols {get('control', 'final_n')}",
         S.RED_FILL, S.RED),
    ]
    for y, h, t, b, fill, edge in boxes:
        S.box(ax, x, y, w, h, t, b, fill=fill, edge=edge, lw=1.1)
    for (y1, *_), (y2, h2, *_) in zip(boxes[:-1], boxes[1:]):
        S.arrow(ax, cx, y1, cx, y2 + h2)

    ex = [("QRS onset before P offset", "qrs_onset_before_p_offset"),
          ("QRS out of range", "range:qrs_ms"),
          ("P out of range", "range:p_ms"),
          ("QT out of range", "range:qt_ms")]
    body = "\n".join(f"{lab}: {get('all', r, 'first')}" for lab, r in ex)
    ex_x, ex_y, ex_w, ex_h = 1.75, 0.78, 1.72, 0.84
    S.box(ax, ex_x, ex_y, ex_w, ex_h, f"Excluded by QC: {get('all', 'excluded')}",
          body, fill="white", edge=S.MUTED, lw=0.9)
    S.arrow(ax, cx, ex_y + ex_h / 2, ex_x, ex_y + ex_h / 2, lw=0.8)
    return fig


# ===========================================================================
# FIG 5 — Bland-Altman, P duration vs LUDB
# ===========================================================================
def fig05_bland_altman(out):
    b = pd.read_csv(out / "gate2_beats.csv")
    t = b[(b.dataset == "ludb") & (b.half == "test") & (b.reference == "knn") & (b.sinus == 1)]
    t = t.dropna(subset=["p_onset_err_ms", "p_offset_err_ms", "p_dur_ref_ms"]).copy()
    t["err"] = t.p_offset_err_ms - t.p_onset_err_ms          # measured - reference
    t["mean"] = t.p_dur_ref_ms + t.err / 2
    bias, sd = t.err.mean(), t.err.std(ddof=1)
    lo, hi = bias - 1.96 * sd, bias + 1.96 * sd

    fig = S.plt.figure(figsize=(S.COL_1, 2.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.2, 1], wspace=0.05)
    ax, axh = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    ax.axhspan(lo, hi, color=S.RED_FILL, zorder=0)
    ax.scatter(t["mean"], t.err, s=7, color=S.RED, alpha=0.45, linewidths=0, zorder=2)
    ax.axhline(0, color=S.MUTED, lw=0.7, zorder=1)
    ax.axhline(bias, color=S.INK, lw=1.1, zorder=3)
    for v in (lo, hi):
        ax.axhline(v, color=S.RED, lw=0.8, ls="--", zorder=3)
    ax.set_xlabel("Mean of measured and expert P duration (ms)")
    ax.set_ylabel("Measured − expert (ms)")
    clean(ax)

    bins = np.arange(np.floor(t.err.min() / 10) * 10, t.err.max() + 10, 10)
    axh.hist(t.err, bins=bins, orientation="horizontal", color=S.RED_FILL,
             edgecolor=S.RED, lw=0.6)
    axh.set_ylim(ax.get_ylim())
    axh.axhline(bias, color=S.INK, lw=1.1)
    axh.set_yticks([])
    axh.set_xticks([])
    for s in ("top", "right", "left", "bottom"):
        axh.spines[s].set_visible(False)
    for v, lab, va in [(hi, f"+1.96 SD\n{hi:+.0f}", "bottom"), (bias, f"bias\n{bias:+.1f}", "center"),
                       (lo, f"−1.96 SD\n{lo:+.0f}", "top")]:
        axh.text(1.02, v, lab, transform=axh.get_yaxis_transform(), fontsize=7.5,
                 va=va, ha="left", color=S.INK, linespacing=1.1)
    fig.subplots_adjust(left=0.17, right=0.80, bottom=0.17, top=0.97)
    return fig


# ===========================================================================
# FIG 6 — case vs control distributions: the headline finding in one picture
# ===========================================================================
PANELS = [("p_duration_ms", "ms"), ("ptfv1_uv_ms", "µV·ms"), ("p_axis_deg", "°"),
          ("qrs_duration_ms", "ms"), ("qtc_fridericia_ms", "ms"), ("st_j60_v6_mv", "mV")]


def _half_violin(ax, x, vals, side, colour, fill, width=0.36):
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    lo, hi = np.percentile(vals, [1, 99])
    grid = np.linspace(lo, hi, 200)
    dens = gaussian_kde(vals)(grid)
    dens = dens / dens.max() * width
    xs = x + side * dens
    ax.fill_betweenx(grid, x, xs, color=fill, lw=0, zorder=1)
    ax.plot(xs, grid, color=colour, lw=0.9, zorder=2)
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    bx = x + side * 0.05
    ax.plot([bx, bx], [q1, q3], color=colour, lw=3.2, solid_capstyle="butt", zorder=3)
    ax.plot(bx, med, marker="o", ms=4, mfc="white", mec=colour, mew=1.1, zorder=4)
    return lo, hi


def fig06_distributions(out):
    q = pd.read_csv(out / "qc_designA.csv")
    q = q[q.qc_pass]
    g = pd.read_csv(out / "group_comparison_designA.csv").set_index("feature")

    fig, axes = S.plt.subplots(2, 3, figsize=(S.COL_2, 4.3))
    rng = np.random.default_rng(0)
    for ax, (f, unit) in zip(axes.flat, PANELS):
        col, fill = axis_colour(f)
        case = q.loc[q.group == "case", f].dropna()
        ctrl = q.loc[q.group == "control", f].dropna()
        lims = []
        for x, vals, c, fl, side in [(0, ctrl, S.MUTED, S.GREY_FILL, -1),
                                     (1, case, col, fill, +1)]:
            lims += _half_violin(ax, x, vals, side, c, fl)
            jit = x - side * (0.10 + 0.16 * rng.random(len(vals)))
            ax.scatter(jit, vals, s=2.2, color=c, alpha=0.35, linewidths=0, zorder=1)
        lo, hi = min(lims), max(lims)
        span = hi - lo
        ax.set_ylim(lo - 0.05 * span, hi + 0.30 * span)
        ax.set_xlim(-0.65, 1.65)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Control\n{len(ctrl)}", f"LAE\n{len(case)}"])
        ax.set_ylabel(unit)
        r = g.loc[f]
        stat = "p" if r.role == "primary" else "q"
        val = r.p if r.role == "primary" else r.q_bh
        sval = f"{stat} = {val:.2f}" if val >= 0.01 else f"{stat} < 0.001" if val < 0.001 else f"{stat} = {val:.3f}"
        ax.text(0.5, 0.97, f"g = {r.hedges_g:+.2f} [{r.g_lo:+.2f}, {r.g_hi:+.2f}]\n{sval}",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                color=S.INK, linespacing=1.15)
        title = flabel(f) + ("  (primary)" if r.role == "primary" else "")
        ax.set_title(title, fontsize=8.5, loc="left", color=col if col != S.GREY else S.INK,
                     fontweight="bold")
        clean(ax)
        ax.tick_params(axis="x", length=0)
    fig.subplots_adjust(hspace=0.62, wspace=0.42, left=0.07, right=0.99, top=0.94, bottom=0.11)
    return fig


# ===========================================================================
# FIG 7 — effect sizes, every marker
# ===========================================================================
GROUPS = [
    ("Atrial", ["p_duration_ms", "ptfv1_uv_ms", "p_axis_deg"]),
    ("Rate", ["mean_hr_bpm", "rmssd_ms"]),
    ("Ventricular", ["qrs_duration_ms", "qt_ms", "qtc_bazett_ms", "qtc_fridericia_ms",
                     "tpe_ms", "tpe_qt_ratio"]),
    ("ST J+60", ["st_j60_i_mv", "st_j60_avl_mv", "st_j60_v5_mv", "st_j60_v6_mv",
                 "st_j60_ii_mv", "st_j60_iii_mv", "st_j60_avf_mv", "st_j60_avr_mv",
                 "st_j60_v1_mv", "st_j60_v2_mv", "st_j60_v3_mv", "st_j60_v4_mv"]),
]


def fig07_forest(out):
    g = pd.read_csv(out / "group_comparison_designA.csv").set_index("feature")
    fig, ax = S.plt.subplots(figsize=(S.COL_1, 4.7))
    fig.subplots_adjust(left=0.36, right=0.97, top=0.98, bottom=0.19)
    y, ticks, heads, bands = 0, [], [], []
    for name, feats in GROUPS:
        heads.append((y, name))
        y0 = y
        y += 1
        for f in feats:
            r = g.loc[f]
            col, _ = axis_colour(f)
            circ = r.circularity == "circular"
            sig = (r.p < 0.05) if r.role == "primary" else (r.q_bh < 0.05)
            ax.plot([r.g_lo, r.g_hi], [y, y], color=col, lw=1.2, solid_capstyle="round")
            ax.plot(r.hedges_g, y, marker="s" if circ else "o", ms=4.6,
                    mfc=col if sig else "white", mec=col, mew=1.1, zorder=3)
            ticks.append((y, flabel(f)))
            y += 1
        bands.append((y0 - 0.5, y - 0.5))
        y += 0.4
    for i, (a, b) in enumerate(bands):
        if i % 2 == 0:
            ax.axhspan(a, b, color=S.GREY_FILL, zorder=0, lw=0)
    ax.axvline(0, color=S.INK_SOFT, lw=0.8)
    ax.set_yticks([t for t, _ in ticks])
    ax.set_yticklabels([l for _, l in ticks], fontsize=7.5)
    for yy, name in heads:
        ax.text(-0.02, yy, name, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=8, fontweight="bold", color=S.INK)
    ax.invert_yaxis()
    ax.set_xlabel("Hedges' g, LAE − control (95% CI)")
    clean(ax, left=False)
    ax.tick_params(axis="y", length=0)
    h = [Line2D([], [], marker="o", ls="", mfc=S.RED, mec=S.RED, ms=4.6, label="Atrial marker"),
         Line2D([], [], marker="o", ls="", mfc=S.BLUE, mec=S.BLUE, ms=4.6, label="Ventricular marker"),
         Line2D([], [], marker="s", ls="", mfc=S.INK_SOFT, mec=S.INK_SOFT, ms=4.6,
                label="LAE criterion (circular)"),
         Line2D([], [], marker="o", ls="", mfc="white", mec=S.INK_SOFT, ms=4.6,
                label="Not significant")]
    ax.legend(handles=h, fontsize=7.5, frameon=False, loc="upper center",
              bbox_to_anchor=(0.4, -0.085), ncol=2, columnspacing=1.0, handletextpad=0.3)
    return fig


# ===========================================================================
# FIG 8 — ROC: best model highlighted, others as context
# ===========================================================================
def fig08_roc(out):
    o = pd.read_csv(out / "models_oof_predictions.csv")
    m = pd.read_csv(out / "models_metrics.csv").set_index(["set", "model"])
    perm = pd.read_csv(out / "models_permutation_test.csv").set_index("set")
    best = "random_forest"

    fig, axes = S.plt.subplots(1, 2, figsize=(S.COL_2, 3.55),
                               gridspec_kw=dict(width_ratios=[1, 1], wspace=0.35))

    # left: both feature sets, best model only, with the permutation null ceiling
    ax = axes[0]
    for st, col, lab in [("A", S.INK, "A · all markers"), ("B", S.BLUE, "B · no P-wave markers")]:
        d = o[(o.set == st) & (o.model == best)]
        fpr, tpr, _ = roc_curve(d.case, d.prob)
        r = m.loc[(st, best)]
        ax.plot(fpr, tpr, color=col, lw=1.8,
                label=f"{lab}\nAUROC {r.auroc:.3f} ({r.auroc_lo:.3f}–{r.auroc_hi:.3f})")
        if st == "A":
            ax.fill_between(fpr, tpr, fpr, color=S.GREY_FILL, zorder=0)
    ax.plot([0, 1], [0, 1], color=S.MUTED, lw=0.8, ls="--")
    nmax = perm["null_max"].max()
    ax.text(0.97, 0.04, f"200 label permutations: AUROC ≤ {nmax:.2f}, p = {perm.p_value.max():.3f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=S.INK_SOFT)
    ax.set_title("Random forest, two feature sets", fontsize=8.5, loc="left", fontweight="bold")
    ax.legend(fontsize=7.5, frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.36),
              handlelength=1.4, labelspacing=0.7)

    # right: all four models per set, as a dot plot of AUROC with CIs
    ax2 = axes[1]
    ys = np.arange(len(MODEL_LABEL))
    for st, col, off, mk in [("A", S.INK, -0.14, "o"), ("B", S.BLUE, +0.14, "D")]:
        for i, mdl in enumerate(MODEL_LABEL):
            r = m.loc[(st, mdl)]
            ax2.plot([r.auroc_lo, r.auroc_hi], [i + off] * 2, color=col, lw=1.2)
            ax2.plot(r.auroc, i + off, marker=mk, ms=4.6, color=col,
                     mfc=col if mdl == best else "white", mew=1.1)
    ax2.axvline(0.5, color=S.MUTED, lw=0.8, ls="--")
    ax2.text(0.505, -0.55, "chance", fontsize=7.5, color=S.INK_SOFT, va="center")
    ax2.set_yticks(ys)
    ax2.set_yticklabels([MODEL_LABEL[k] for k in MODEL_LABEL])
    ax2.set_ylim(len(MODEL_LABEL) - 0.5, -0.8)
    ax2.set_xlim(0.45, 0.9)
    ax2.set_xlabel("AUROC (95% CI), nested CV")
    ax2.set_title("All models", fontsize=8.5, loc="left", fontweight="bold")
    ax2.legend(handles=[Line2D([], [], marker="o", color=S.INK, ms=4.6, label="Set A"),
                        Line2D([], [], marker="D", color=S.BLUE, ms=4.6, label="Set B")],
               fontsize=7.5, frameon=False, loc="lower left")
    ax2.tick_params(axis="y", length=0)

    for a in axes:
        clean(a)
    axes[0].set_xlim(0, 1); axes[0].set_ylim(0, 1); axes[0].set_aspect("equal")
    axes[0].set_xlabel("1 − specificity")
    axes[0].set_ylabel("Sensitivity")
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.14, top=0.92)
    return fig


# ===========================================================================
# FIG 9 — permutation importance, paired A vs B
# ===========================================================================
def fig09_importance(out, k=10):
    imp = pd.read_csv(out / "models_importance.csv")
    imp = imp[imp.model == "random_forest"]
    fig, axes = S.plt.subplots(1, 2, figsize=(S.COL_2, 2.8),
                               gridspec_kw=dict(wspace=0.85))
    for ax, st, title in zip(axes, ["A", "B"], ["A · all markers", "B · no P-wave markers"]):
        d = imp[imp.set == st].sort_values("importance_mean", ascending=False).head(k)[::-1]
        cols = [axis_colour(f) for f in d.feature]
        ax.barh(range(len(d)), d.importance_mean, height=0.62,
                color=[c[1] for c in cols], edgecolor=[c[0] for c in cols], lw=0.9)
        ax.errorbar(d.importance_mean, range(len(d)), xerr=d.importance_sd_across_folds,
                    fmt="none", ecolor=S.INK_SOFT, elinewidth=0.7, capsize=1.6)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels([flabel(f) for f in d.feature])
        ax.set_xlabel("Drop in AUROC when shuffled")
        ax.xaxis.set_major_locator(S.plt.MaxNLocator(4))
        ax.set_title(f"{title} (random forest)", fontsize=8.5, loc="left", fontweight="bold")
        ax.tick_params(axis="y", length=0)
        clean(ax)
    axes[1].legend(handles=[S.plt.Rectangle((0, 0), 1, 1, fc=S.RED_FILL, ec=S.RED, label="Atrial"),
                            S.plt.Rectangle((0, 0), 1, 1, fc=S.BLUE_FILL, ec=S.BLUE, label="Ventricular"),
                            S.plt.Rectangle((0, 0), 1, 1, fc=S.GREY_FILL, ec=S.GREY, label="Rate / age")],
                   fontsize=7.5, frameon=False, loc="lower right")
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.17, top=0.9)
    return fig


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--fig-dir", default="fig")
    a = ap.parse_args()
    for name, fn in [("fig04_cohort_flow", fig04_flow),
                     ("fig05_bland_altman", fig05_bland_altman),
                     ("fig06_distributions", fig06_distributions),
                     ("fig07_effect_sizes", fig07_forest),
                     ("fig08_roc", fig08_roc),
                     ("fig09_importance", fig09_importance)]:
        S.save(fn(a.out_dir), name, a.fig_dir)
