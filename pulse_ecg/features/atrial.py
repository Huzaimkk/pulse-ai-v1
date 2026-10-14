"""
Atrial family — the atrial-cardiopathy axis.

P duration is global: multilead P onset to P offset on the median beat (Rule 3:
all leads present, never one chosen lead). It is the PRIMARY ENDPOINT
(config.yaml analysis.primary). PTFV1 is read in V1 on V1's own single-lead P
marks (delineation.per_lead), not the pooled ones. P axis, IAB morphology and the
MVP voltage term are read on the median beat inside the global P window; the
median beat is baseline-referenced to the isoelectric run before QRS onset
(delineate.py).

Every normal range and threshold is read from config.yaml `atrial:`.
"""

import math

import numpy as np

from .. import config
from ..registry import NotApplicable, feature

A = config.load()["atrial"]
IAB = A["iab"]
MVP = A["mvp"]
MIN_S = A["min_duration_s"]
MIN_MV = A["min_deflection_mv"]
INFERIOR = list(IAB["inferior_leads"])


def _thr(entry):
    return tuple(entry["abnormal"]) if "abnormal" in entry else None


# ------------------------------------------------------------- helpers ---
def _idx(ctx, samples: float) -> int:
    """Median-beat sample index, or NotApplicable if outside the beat window."""
    if samples is None or math.isnan(samples):
        raise NotApplicable("P boundary not delineated")
    i = int(round(samples))
    if not 0 <= i < len(ctx.median_beat):
        raise NotApplicable("point falls outside the median-beat window")
    return i


def _p_window(ctx, lead: str, onset: float, offset: float) -> np.ndarray:
    lo, hi = _idx(ctx, onset), _idx(ctx, offset)
    if hi <= lo:
        raise NotApplicable(f"P offset does not follow P onset in {lead}")
    return np.asarray(ctx.median_beat[lo:hi + 1, ctx.rec.leads.index(lead)], dtype=float)


def _global_p(ctx, lead: str) -> np.ndarray:
    """Lead `lead` of the median beat inside the multilead P window."""
    return _p_window(ctx, lead, ctx.points["p_onset"], ctx.points["p_offset"])


def terminal_negative(x, fs: float, min_mv: float = MIN_MV) -> tuple[float, float]:
    """(duration ms, depth µV) of the terminal negative part of a P wave: the run
    of negative samples around the P trough, bounded by the nearest non-negative
    sample on each side. (0, 0) if the trough does not reach -min_mv (an upright
    P) or a positive deflection >= min_mv follows it (not terminal). Sub-min_mv
    baseline ripple after the wave neither ends nor disqualifies the component."""
    x = np.asarray(x, dtype=float)
    i = int(np.argmin(x))
    if x[i] > -min_mv:
        return 0.0, 0.0
    before = np.flatnonzero(x[:i] >= 0)
    after = np.flatnonzero(x[i:] >= 0)
    lo = before[-1] + 1 if len(before) else 0
    hi = i + after[0] if len(after) else len(x)
    if hi < len(x) and x[hi:].max() >= min_mv:
        return 0.0, 0.0
    return (hi - lo) * 1000.0 / fs, -float(x[i]) * 1000.0


def biphasic_pos_neg(x, min_mv: float = MIN_MV) -> bool:
    """+/- P: a positive deflection >= min_mv followed by a negative one <= -min_mv."""
    x = np.asarray(x, dtype=float)
    i_max, i_min = int(np.argmax(x)), int(np.argmin(x))
    return bool(x[i_max] >= min_mv and x[i_min] <= -min_mv and i_max < i_min)


def p_axis_from_areas(area_i: float, area_avf: float) -> float:
    """Frontal-plane axis in degrees. Hexaxial geometry (Einthoven/Goldberger):
    I = E cos(theta), aVF = E (sqrt(3)/2) sin(theta)."""
    return math.degrees(math.atan2(2.0 * area_avf / math.sqrt(3.0), area_i))


def classify_iab(p_ms: float, n_biphasic_inferior: int) -> int:
    """0 none, 1 partial (P >= iab.p_ms), 2 advanced (and +/- P in
    >= iab.advanced_min_leads inferior leads)."""
    if p_ms < IAB["p_ms"]:
        return 0
    return 2 if n_biphasic_inferior >= IAB["advanced_min_leads"] else 1


def mvp_points(iab: int, p_lead_i_mv: float, p_ms: float) -> int:
    """MVP ECG risk score: Morphology (IAB class 0-2) + Voltage (P amplitude in
    lead I) + P duration, each 0-2 points."""
    v_lo, v_hi = MVP["voltage_mv"]
    d_lo, d_hi = MVP["duration_ms"]
    v = 0 if p_lead_i_mv > v_hi else (1 if p_lead_i_mv >= v_lo else 2)
    d = 0 if p_ms < d_lo else (1 if p_ms <= d_hi else 2)
    return int(iab) + v + d


# ------------------------------------------------------------ P duration ---
@feature(name="p_duration_ms", family="atrial", units="ms", leads=["vm"],
         requires=["p_onset", "p_offset"], min_duration_s=MIN_S,
         normal=A["p_duration_ms"]["normal"], abnormal=_thr(A["p_duration_ms"]),
         cite="martinezselles2020bayes",
         desc="Global P-wave duration, multilead P onset to P offset (12-lead vector "
              "magnitude). Primary endpoint, analysed as a continuous variable")
def p_duration_ms(ctx):
    """PRIMARY ENDPOINT. Global P-wave duration in ms: multilead P onset to P
    offset over every lead present. The paper's main result is the continuous
    value; the 120 ms threshold is for Table III and descriptive use only."""
    return (ctx.points["p_offset"] - ctx.points["p_onset"]) * 1000.0 / ctx.fs


# ----------------------------------------------------------------- PTFV1 ---
@feature(name="ptfv1_uv_ms", family="atrial", units="µV·ms", leads=["V1"],
         requires=["median_beat"], min_duration_s=MIN_S,
         normal=A["ptfv1_uv_ms"]["normal"], abnormal=_thr(A["ptfv1_uv_ms"]),
         cite="morris1964ptfv1",
         desc="P-terminal force in V1 (Morris index): -(depth in µV) x (duration in ms) "
              "of the terminal negative P component, on V1's own P boundaries")
def ptfv1_uv_ms(ctx):
    """Morris P-terminal force in V1, µV·ms, on V1's own single-lead P marks.
    Negative for a biphasic V1 P, 0 for an upright one."""
    if "V1" not in ctx.rec.leads:
        raise NotApplicable("requires lead V1, not present (single-lead device path)")
    marks = (ctx.extra.get("per_lead") or {}).get("V1") or {}
    on, off = marks.get("p_onset", math.nan), marks.get("p_offset", math.nan)
    if any(v is None or math.isnan(v) for v in (on, off)):
        raise NotApplicable("no single-lead P boundaries in V1")
    dur_ms, depth_uv = terminal_negative(_p_window(ctx, "V1", on, off), ctx.fs)
    return -depth_uv * dur_ms


# ---------------------------------------------------------------- P axis ---
@feature(name="p_axis_deg", family="atrial", units="deg", leads=["I", "aVF"],
         requires=["p_onset", "p_offset", "median_beat"], min_duration_s=MIN_S,
         normal=A["p_axis_deg"]["normal"], abnormal=None, cite="maheshwari2017paxis",
         desc="Frontal-plane P axis from the P-wave areas in I and aVF (global P window)")
def p_axis_deg(ctx):
    return p_axis_from_areas(float(np.sum(_global_p(ctx, "I"))),
                             float(np.sum(_global_p(ctx, "aVF"))))


# ------------------------------------------------------------------- IAB ---
@feature(name="iab_class", family="atrial", units="class", leads=INFERIOR,
         requires=["p_onset", "p_offset", "median_beat"], min_duration_s=MIN_S,
         normal=None, abnormal=_thr(IAB), cite="bayesdeluna2012iab",
         desc="Interatrial block: 0 none, 1 partial (P >= 120 ms), 2 advanced (and "
              "+/- P in an inferior lead). Secondary, descriptive only")
def iab_class(ctx):
    """SECONDARY / DESCRIPTIVE ONLY — never an endpoint. 0 = none, 1 = partial
    (P >= 120 ms), 2 = advanced (P >= 120 ms and biphasic +/- P in II, III or aVF).

    The binary call is not trustworthy: on LUDB sinus records 66% measure
    >= 120 ms under any pooling rule, against a true IAB prevalence of 10-20%.
    The continuous p_duration_ms is the endpoint."""
    n = sum(biphasic_pos_neg(_global_p(ctx, l)) for l in INFERIOR)
    return classify_iab(p_duration_ms(ctx), n)


# ------------------------------------------------------------------- MVP ---
@feature(name="mvp_score", family="atrial", units="points", leads=["I", *INFERIOR],
         requires=["p_onset", "p_offset", "median_beat"], min_duration_s=MIN_S,
         normal=MVP["normal"], abnormal=_thr(MVP), cite="alexander2019mvp",
         desc="MVP ECG risk score (0-6): P morphology (IAB class) + P voltage in "
              "lead I + P duration")
def mvp_score(ctx):
    """Morphology-Voltage-P-wave duration score, 0-6. The morphology term is
    iab_class, so it inherits that feature's descriptive-only caveat."""
    v = float(np.max(np.abs(_global_p(ctx, "I"))))
    return mvp_points(iab_class(ctx), v, p_duration_ms(ctx))
