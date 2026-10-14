"""
Ventricular family — the arrhythmia axis.

Global intervals (QRS, QT, QTc, Tpe) come from the multilead delineation points on
the median beat (Rule 3: all leads present, never one chosen lead). fQRS and the
ST level are per lead, on the lead's own median beat, baseline-referenced to the
isoelectric run before QRS onset (delineate.py).

Every normal range and threshold is read from config.yaml `ventricular:`.
"""

import math

import numpy as np
from scipy import signal as ss

from .. import config
from ..delineate import median_beat
from ..preprocess import filter_signal
from ..recording import STANDARD_LEADS
from ..registry import NotApplicable, feature
from .rate import rr_ms

V = config.load()["ventricular"]
FQRS = V["fqrs"]
QTC = V["qtc_ms"]


def _thr(entry):
    return tuple(entry["abnormal"]) if "abnormal" in entry else None


# ------------------------------------------------------------- helpers ---
def _ms(ctx, a, b) -> float:
    return (ctx.points[b] - ctx.points[a]) * 1000.0 / ctx.fs


def _idx(ctx, samples: float) -> int:
    """Median-beat sample index, or NotApplicable if outside the beat window."""
    i = int(round(samples))
    if not 0 <= i < len(ctx.median_beat):
        raise NotApplicable("point falls outside the median-beat window")
    return i


def _rr_s(ctx) -> float:
    rr = rr_ms(ctx.r_peaks, ctx.fs)
    if len(rr) < 1:
        raise NotApplicable("fewer than 2 detected beats")
    return float(np.mean(rr)) / 1000.0


def _t_peak(ctx) -> int:
    """Spatial T peak: maximum of the vector magnitude between J + st_after_j_ms
    and T offset."""
    lo = _idx(ctx, ctx.points["qrs_offset"] + V["st_after_j_ms"] * ctx.fs / 1000)
    hi = _idx(ctx, ctx.points["t_offset"])
    if hi <= lo:
        raise NotApplicable("T offset precedes J + st_after_j_ms")
    return lo + int(np.argmax(ctx.vm[lo:hi + 1]))


def fragmented(x, min_prominence_mv: float | None = None) -> bool:
    """Das et al. 2006 fQRS morphology in one lead's QRS: an additional R wave
    (R') or a notch in the nadir of the S (or Q) wave.

    Deflections are local extrema with prominence >= min_prominence_mv.
      * >= 2 positive peaks (rSR', notched R): fragmented.
      * two troughs separated by a peak smaller in magnitude than both (a notched
        S or Q nadir): fragmented. A qRs has two troughs too, but separated by
        the R wave, which is larger than both, so it is not."""
    p = FQRS["min_prominence_mv"] if min_prominence_mv is None else min_prominence_mv
    x = np.asarray(x, dtype=float)
    peaks, _ = ss.find_peaks(x, prominence=p)
    if len(peaks) >= 2:
        return True
    troughs, _ = ss.find_peaks(-x, prominence=p)
    for a, b in zip(troughs[:-1], troughs[1:]):
        if abs(np.max(x[a:b + 1])) < min(abs(x[a]), abs(x[b])):
            return True
    return False


# ------------------------------------------------------ global intervals ---
@feature(name="qrs_duration_ms", family="ventricular", units="ms", leads=["vm"],
         requires=["qrs_onset", "qrs_offset"], min_duration_s=0,
         normal=V["qrs_duration_ms"]["normal"], abnormal=_thr(V["qrs_duration_ms"]),
         cite="surawicz2009aha3",
         desc="Global QRS duration, multilead QRS onset to QRS offset")
def qrs_duration_ms(ctx):
    return _ms(ctx, "qrs_onset", "qrs_offset")


@feature(name="qt_ms", family="ventricular", units="ms", leads=["vm"],
         requires=["qrs_onset", "t_offset"], min_duration_s=0,
         normal=V["qt_ms"]["normal"], abnormal=_thr(V["qt_ms"]), cite="rautaharju2009aha4",
         desc="Global QT interval, multilead QRS onset to T offset")
def qt_ms(ctx):
    return _ms(ctx, "qrs_onset", "t_offset")


@feature(name="qtc_bazett_ms", family="ventricular", units="ms", leads=["vm"],
         requires=["qrs_onset", "t_offset", "r_peaks"], min_duration_s=0,
         normal=QTC["normal"], abnormal=_thr(QTC),
         cite=["bazett1920qt", "rautaharju2009aha4"],
         desc="QT corrected by Bazett, QT / sqrt(mean RR in s)")
def qtc_bazett_ms(ctx):
    return qt_ms(ctx) / math.sqrt(_rr_s(ctx))


@feature(name="qtc_fridericia_ms", family="ventricular", units="ms", leads=["vm"],
         requires=["qrs_onset", "t_offset", "r_peaks"], min_duration_s=0,
         normal=QTC["normal"], abnormal=_thr(QTC),
         cite=["fridericia1920qt", "rautaharju2009aha4"],
         desc="QT corrected by Fridericia, QT / cbrt(mean RR in s)")
def qtc_fridericia_ms(ctx):
    return qt_ms(ctx) / _rr_s(ctx) ** (1.0 / 3.0)


@feature(name="tpe_ms", family="ventricular", units="ms", leads=["vm"],
         requires=["qrs_offset", "t_offset", "vm"], min_duration_s=0,
         normal=None, abnormal=None, cite="gupta2008tpe",
         desc="T peak to T end; spatial T peak (vector-magnitude maximum) to multilead T offset")
def tpe_ms(ctx):
    return (ctx.points["t_offset"] - _t_peak(ctx)) * 1000.0 / ctx.fs


@feature(name="tpe_qt_ratio", family="ventricular", units="", leads=["vm"],
         requires=["qrs_onset", "qrs_offset", "t_offset", "vm"], min_duration_s=0,
         normal=V["tpe_qt_ratio"]["normal"], abnormal=_thr(V["tpe_qt_ratio"]),
         cite="gupta2008tpe", desc="Tpeak-Tend / QT")
def tpe_qt_ratio(ctx):
    return tpe_ms(ctx) / qt_ms(ctx)


# ------------------------------------------------------------------ fQRS ---
_FQRS_LEADS = [l for group in FQRS["territories"].values() for l in group]


def _wideband_beat(ctx) -> np.ndarray:
    """Median beat of the unfiltered recording, band-passed to fqrs.band_hz instead
    of the 40 Hz analysis band (which smooths notches away). Same beat window and
    R-peaks as ctx.median_beat, so delineation points index it directly; baseline is
    the isoelectric run before QRS onset, as in delineate.py."""
    lo_hz, hi_hz = FQRS["band_hz"]
    hi_hz = min(hi_hz, 0.45 * ctx.fs)
    cfg = {**ctx.cfg, "filter": {**ctx.cfg["filter"], "band": [lo_hz, hi_hz]}}
    x = filter_signal(ctx.raw.signal, ctx.fs, cfg)
    med, r0, _, _ = median_beat(x, ctx.r_peaks, ctx.fs, ctx.cfg)
    if med is None or r0 != ctx.r0:
        raise NotApplicable("no wide-band median beat")
    run = max(2, int(round(ctx.cfg["bounds"]["min_run_ms"] * ctx.fs / 1000)))
    q = _idx(ctx, ctx.points["qrs_onset"])
    return med - med[max(0, q - run):max(q, 1)].mean(axis=0)


@feature(name="fqrs_present", family="ventricular", units="", leads=_FQRS_LEADS,
         requires=["qrs_onset", "qrs_offset", "median_beat"], min_duration_s=0,
         normal=None, abnormal=("==", 1), cite="das2006fqrs",
         desc="Fragmented QRS (Das criteria): RSR' or notched R/S/Q nadir in >= 2 leads "
              "of one coronary territory; narrow QRS only")
def fqrs_present(ctx):
    if qrs_duration_ms(ctx) >= FQRS["max_qrs_ms"]:
        raise NotApplicable(f"QRS >= {FQRS['max_qrs_ms']} ms; Das narrow-QRS criteria do not apply")
    lo, hi = _idx(ctx, ctx.points["qrs_onset"]), _idx(ctx, ctx.points["qrs_offset"])
    beat = _wideband_beat(ctx)
    frag = {l: fragmented(beat[lo:hi + 1, ctx.rec.leads.index(l)]) for l in _FQRS_LEADS}
    hit = any(sum(frag[l] for l in group) >= FQRS["min_leads"]
              for group in FQRS["territories"].values())
    return int(hit)


# ----------------------------------------------------------- ST at J+60 ---
def _st_feature(lead):
    name = f"st_j60_{lead.lower()}_mv"

    @feature(name=name, family="ventricular", units="mV", leads=[lead],
             requires=["qrs_offset", "median_beat"], min_duration_s=0,
             normal=None, abnormal=None, cite="rautaharju2009aha4",
             desc=f"ST level in {lead} at J + {V['st_after_j_ms']} ms, relative to the "
                  f"isoelectric run before QRS onset")
    def st(ctx):
        i = _idx(ctx, ctx.points["qrs_offset"] + V["st_after_j_ms"] * ctx.fs / 1000)
        return float(ctx.median_beat[i, ctx.rec.leads.index(lead)])
    return st


for _lead in STANDARD_LEADS:
    _st_feature(_lead)
