"""
make_figures.py — every figure that needs no data
==================================================
Generates Figs. 1-3 for the PULSE AI paper. None of these depend on results,
so they can be finalised now and will not change when the numbers arrive.

    python make_figures.py

Output: fig/fig01_system.pdf, fig/fig02_pipeline.pdf, fig/fig03_chain.pdf
(plus PNGs for quick viewing; the PDFs are what go into LaTeX).

All coordinates are inches (see figstyle.canvas). S.save() refuses to write a
figure whose text is below 7.5 pt, overlaps, or escapes its box.
"""

import figstyle as S

S.apply()

W = S.COL_2
M = 0.04          # outer margin, inches


# ===========================================================================
# FIG 1 — Wearable prototype, layered block diagram
# ===========================================================================
def fig01_system():
    from diagram_parts import (conn, icon_cap, icon_chip, icon_electrodes, icon_optical,
                               icon_screen, item, layer)
    H = 2.95
    fig, ax = S.canvas(W, H)
    x0, lw_ = M, W - 2 * M
    LB = 0.95                         # label column width
    cx0 = x0 + LB + 0.12              # content starts here

    # bottom: sensing layer
    layer(ax, x0, 0.04, lw_, 0.82, "Sensing\nlayer", LB)
    sw, sg = 1.80, 0.17
    sens = [("AD8232", "single-lead ECG", icon_electrodes),
            ("MAX30102", "PPG, SpO$_2$", icon_optical),
            ("OpenBCI cap", "EEG (later papers)", icon_cap)]
    sx = [cx0 + i * (sw + sg) for i in range(3)]
    for x, (t, b, ic) in zip(sx, sens):
        item(ax, x, 0.16, sw, 0.58, t, b, icon=ic)

    # middle: edge layer
    layer(ax, x0, 1.06, lw_, 0.82, "Edge\nlayer", LB)
    ex, ew = cx0 + 1.30, 3.00
    item(ax, ex, 1.18, ew, 0.58, "ESP32 edge node", "12-bit ADC sampling, time stamps,\nMQTT publish over Wi-Fi",
         icon=icon_chip)

    # top: processing and output
    layer(ax, x0, 2.08, lw_, 0.83, "Processing\n& output", LB)
    px, pw = cx0, 2.70
    item(ax, px, 2.20, pw, 0.59, "Raspberry Pi 3B+", "filtering, median beat, marker\nextraction (same code as Fig. 2)",
         icon=icon_chip, edge=S.INK, lw=1.1)
    ox, ow = px + pw + 0.55, 2.05
    item(ax, ox, 2.20, ow, 0.59, "OLED / dashboard", "markers and trends\nfor the user", icon=icon_screen)

    for i, x in enumerate(sx):
        conn(ax, x + sw / 2, 0.74, ex + 0.6 + i * 0.9, 1.18, f"1.{i + 1}" if i != 1 else "1.2 sample",
             side="right")
    conn(ax, ex + ew / 2, 1.76, px + pw / 2 + 0.3, 2.20, "2. MQTT", side="right")
    conn(ax, px + pw, 2.495, ox, 2.495, "3. show")
    return fig


# ===========================================================================
# FIG 2 — Analysis pipeline, layered block diagram
# ===========================================================================
def fig02_pipeline():
    from diagram_parts import conn, icon_db, icon_ecg, item, layer
    H = 4.81
    fig, ax = S.canvas(W, H)
    x0, lw_ = M, W - 2 * M
    LB = 0.95
    cx0 = x0 + LB + 0.12
    span = W - M - 0.10 - cx0

    # ---- layer 1: data
    y1 = 3.92
    layer(ax, x0, y1, lw_, 0.85, "Data", LB)
    dw, dg = (span - 2 * 0.18) / 3, 0.18
    data = [("PTB-XL v1.0.3", "21,799 ECGs, 12-lead,\n500 Hz, 10 s"),
            ("LUDB", "expert wave boundaries,\n12-lead (tuning / held-out)"),
            ("QTDB", "external check,\n2-lead, 250 Hz")]
    dx = [cx0 + i * (dw + dg) for i in range(3)]
    for x, (t, b) in zip(dx, data):
        item(ax, x, y1 + 0.12, dw, 0.61, t, b, icon=icon_db)

    # ---- layer 2: signal processing
    y2 = 2.74
    layer(ax, x0, y2, lw_, 0.98, "Signal\nprocessing", LB)
    steps = [("2.1 Filtering", "50 Hz notch\n0.5–40 Hz, zero phase"),
             ("2.2 Median beat", "R-aligned beats,\nper-sample median"),
             ("2.3 Delineation", "wavelet, per lead\n(Martínez)"),
             ("2.4 Multilead rule", "k-of-N agreement\n(Laguna)")]
    gap = 0.24
    sw = (span - 3 * gap) / 4
    sx = [cx0 + i * (sw + gap) for i in range(4)]
    for i, (x, (t, b)) in enumerate(zip(sx, steps)):
        item(ax, x, y2 + 0.13, sw, 0.62, t, b,
             edge=S.INK if i >= 2 else S.MUTED, lw=1.1 if i >= 2 else 0.9)
        if i < 3:
            conn(ax, x + sw, y2 + 0.44, x + sw + gap, y2 + 0.44, color=S.MUTED)
    icon_ecg(ax, sx[1] + sw - 0.5, y2 + 0.80, 0.42, 0.13)

    conn(ax, dx[0] + dw / 2, y1 + 0.12, sx[0] + sw / 2, y2 + 0.75, "1.1 records500", side="right")
    conn(ax, dx[1] + dw / 2, y1 + 0.12, sx[2] + sw / 2 - 0.25, y2 + 0.75, "1.2 tune / validate",
         dashed=True, color=S.INK_SOFT, side="left")
    conn(ax, dx[2] + dw / 2, y1 + 0.12, sx[2] + sw / 2 + 0.25, y2 + 0.75, "1.3 external",
         dashed=True, color=S.INK_SOFT, side="right")

    # ---- layer 3: markers
    y3 = 1.24
    layer(ax, x0, y3, lw_, 1.30, "Marker\nextraction", LB)
    mg = 0.18
    mw = (span - 2 * mg) / 3
    mx = [cx0 + i * (mw + mg) for i in range(3)]
    marks = [("Atrial", "P duration, PTFV1 (V1)\nP axis (I, aVF), IAB, MVP", S.RED_FILL, S.RED),
             ("Ventricular", "QRS, QT, QTc (B, F)\nTpeak–Tend, Tpe/QT, fQRS, ST", S.BLUE_FILL, S.BLUE),
             ("Rate", "mean HR, RMSSD\n(no LF/HF or SDNN at 10 s)", S.GREY_FILL, S.MUTED)]
    for x, (t, b, fl, ed) in zip(mx, marks):
        S.box(ax, x, y3 + 0.13, mw, 0.82, t, b, fill=fl, edge=ed, lw=1.1)
    for i, x in enumerate(mx):
        conn(ax, sx[3] + sw / 2 if i == 2 else sx[2] + sw, y2 + 0.13, x + mw / 2, y3 + 0.95,
             f"3.{i + 1}", side="left" if i == 0 else "right")

    # ---- layer 4: analysis
    y4 = 0.04
    layer(ax, x0, y4, lw_, 0.98, "Analysis", LB)
    an = [("4.1 Quality gate", "QRS onset after P offset,\nranges; blind to group"),
          ("4.2 Group comparison", "pre-specified, OLS\n+ age, sex; BH for 23"),
          ("4.3 Classifiers", "4 models, nested CV,\ngrouped by patient"),
          ("Outputs", "effect size per marker,\nAUROC with 95% CI")]
    aw = (span - 3 * gap) / 4
    ax_ = [cx0 + i * (aw + gap) for i in range(4)]
    for i, (x, (t, b)) in enumerate(zip(ax_, an)):
        item(ax, x, y4 + 0.15, aw, 0.66, t, b, edge=S.INK if i == 3 else S.MUTED,
             lw=1.1 if i == 3 else 0.9)
    conn(ax, ax_[0] + aw, y4 + 0.48, ax_[1], y4 + 0.48, color=S.MUTED)
    conn(ax, ax_[1] + aw, y4 + 0.48, ax_[2], y4 + 0.48, color=S.MUTED)
    conn(ax, ax_[2] + aw, y4 + 0.48, ax_[3], y4 + 0.48, color=S.MUTED)
    conn(ax, mx[0] + mw / 2, y3 + 0.13, ax_[0] + aw / 2, y4 + 0.81, "4. markers, one row per record",
         side="right")
    return fig


# ===========================================================================
# FIG 3 — Clinical rationale
# ===========================================================================
def fig03_chain():
    H = 1.72
    fig, ax = S.canvas(W, H)

    y, h, w, gap = 0.52, 1.14, 1.92, 0.66
    stages = [
        ("Stage 1 — substrate", "Atrial cardiopathy",
         "P duration $\\geq$ 120 ms\nadvanced interatrial block\nPTFV1 < −4000 µV·ms\nabnormal P axis",
         S.RED_FILL, S.RED),
        ("Stage 2 — arrhythmia", "Atrial fibrillation",
         "irregular RR intervals\nabsent P waves\noften silent\nmany cases paroxysmal",
         S.GREY_FILL, S.MUTED),
        ("Stage 3 — outcome", "Cardioembolic stroke",
         # ARCADIA (kamel2024arcadia): apixaban in atrial cardiopathy without
         # AF gave no benefit. Never imply anticoagulation helps Stage 1.
         "anticoagulation reduces stroke\nrisk in documented AF; benefit\n"
         "in atrial cardiopathy without\nAF is unproven (ARCADIA)",
         S.GREY_FILL, S.MUTED),
    ]
    xs = [M + i * (w + gap) for i in range(3)]
    for x, (eyebrow, title, body, fill, edge) in zip(xs, stages):
        patch = S.box(ax, x, y, w, h, "", None, fill=fill, edge=edge, lw=1.1)
        tx, top = x + 0.09, y + h
        S._text_in(ax, patch, tx, top, -6, eyebrow, ha="left", va="top",
                   fontsize=S.MIN_PT, color=S.MUTED, fontweight="bold")
        S._text_in(ax, patch, tx, top, -17, title, ha="left", va="top",
                   fontsize=8.6, fontweight="bold", color=S.INK)
        S._text_in(ax, patch, tx, top, -33, body, ha="left", va="top",
                   fontsize=S.BODY_PT, color=S.INK_SOFT, linespacing=1.4)

    for i, (top, bot) in enumerate([("HR 3.31", "incident AF"), ("≈ 5×", "stroke risk")]):
        x1, x2 = xs[i] + w, xs[i + 1]
        mid = (x1 + x2) / 2
        S.arrow(ax, x1 + 0.06, y + h / 2, x2 - 0.06, y + h / 2, lw=1.1)
        ax.text(mid, y + h / 2 + 0.08, top, ha="center", va="bottom",
                fontsize=S.MIN_PT, fontweight="bold", color=S.INK)
        ax.text(mid, y + h / 2 - 0.08, bot, ha="center", va="top",
                fontsize=S.MIN_PT, color=S.MUTED)

    ax.text(M, 0.40,
            "Stage 1 is detectable on a resting ECG before any arrhythmia has occurred.\n"
            # Each number cites its own source:
            #   HR 3.31 -> martinezselles2020bayes, ≈5× -> wolf1991framingham
            "HR 3.31 (95% CI 1.87–5.86): advanced IAB → incident AF, BAYES registry "
            "(Martínez-Sellés 2020).\n"
            "≈5×: stroke risk with AF, Framingham Study (Wolf 1991). "
            "This work measures Stage 1 only.",
            fontsize=S.BODY_PT, color=S.INK_SOFT, va="top", linespacing=1.3)
    return fig


FIGURES = {
    "fig01_system": fig01_system,
    "fig02_pipeline": fig02_pipeline,
    "fig03_chain": fig03_chain,
}

# ===========================================================================
if __name__ == "__main__":
    print("Generating figures that need no data:\n")
    for name, make in FIGURES.items():
        S.save(make(), name)
    print("\nDone. Put the PDFs in your Overleaf project.")
