"""
Rate family — duration-gated.

On every record: n_beats, mean_hr_bpm, rr_intervals_ms, rmssd_ms (the only rate
indices valid on 10 s, Rule 5). Declared with min_duration_s = hrv.long_window_s
(300 s): sdnn_ms, pnn50_pct, lf_ms2, hf_ms2, lf_hf_ratio, sampen, dfa_alpha1. They
stay NaN ("requires >= 300 s") on PTB-XL and fire on long device recordings.

The helpers take an RR series in ms, so they are tested on known series
independently of delineation.
"""

import math

import numpy as np
from scipy import interpolate, signal as ss

from .. import config
from ..registry import NotApplicable, feature

LONG_S = config.load()["hrv"]["long_window_s"]


# ------------------------------------------------------------- helpers ---
def rr_ms(r_peaks, fs) -> np.ndarray:
    return np.diff(np.asarray(r_peaks, dtype=float)) * 1000.0 / fs


def rmssd(rr) -> float:
    return float(np.sqrt(np.mean(np.diff(rr) ** 2)))


def sdnn(rr) -> float:
    return float(np.std(rr, ddof=1))


def pnn50(rr, threshold_ms: float | None = None) -> float:
    thr = config.load()["hrv"]["pnn_ms"] if threshold_ms is None else threshold_ms
    return float(100.0 * np.mean(np.abs(np.diff(rr)) > thr))


def _require_variability(rr, resolution_ms: float = 0.0):
    """Variability below one sample period is quantisation, not physiology."""
    if len(rr) < 3 or np.std(rr) <= resolution_ms:
        raise NotApplicable("RR series has no variability above the sampling resolution")


def lf_hf(rr, resolution_ms: float = 0.0) -> tuple[float, float, float]:
    """(LF ms², HF ms², LF/HF): Welch spectrum of the evenly resampled, linearly
    detrended tachogram, integrated over the Task Force bands."""
    rr = np.asarray(rr, dtype=float)
    _require_variability(rr, resolution_ms)
    h = config.load()["hrv"]
    t = np.cumsum(rr) / 1000.0
    fs = h["resample_hz"]
    grid = np.arange(t[0], t[-1], 1.0 / fs)
    x = ss.detrend(interpolate.CubicSpline(t, rr)(grid))
    f, pxx = ss.welch(x, fs=fs, nperseg=min(len(x), int(h["welch_seg_s"] * fs)))

    def band(lo, hi):
        m = (f >= lo) & (f < hi)
        return float(np.trapezoid(pxx[m], f[m]))
    lf, hf = band(*h["lf_band"]), band(*h["hf_band"])
    if hf <= 0:
        raise NotApplicable("no power in the HF band")
    return lf, hf, lf / hf


def sampen(x, m: int | None = None, r_frac: float | None = None,
           resolution_ms: float = 0.0) -> float:
    """Sample entropy (Richman & Moorman 2000): -ln(A/B); Chebyshev distance <= r,
    r = r_frac x SD, self-matches excluded, the same N - m templates for both lengths."""
    c = config.load()["hrv"]["sampen"]
    m = c["m"] if m is None else m
    x = np.asarray(x, dtype=float)
    _require_variability(x, resolution_ms)
    r = (c["r_frac"] if r_frac is None else r_frac) * float(np.std(x))
    emb = np.lib.stride_tricks.sliding_window_view(x, m + 1)
    a = b = 0
    for i in range(len(emb) - 1):
        close = np.max(np.abs(emb[i + 1:, :m] - emb[i, :m]), axis=1) <= r
        b += int(close.sum())
        a += int((close & (np.abs(emb[i + 1:, m] - emb[i, m]) <= r)).sum())
    if a == 0 or b == 0:
        raise NotApplicable("no template matches; series too short for SampEn")
    return float(-math.log(a / b))


def dfa_alpha1(x, beats=None, resolution_ms: float = 0.0) -> float:
    """Short-term DFA exponent (Peng et al. 1995): slope of log F(n) on log n over
    box sizes n in [lo, hi] beats, linear detrending in non-overlapping boxes."""
    lo, hi = config.load()["hrv"]["dfa_alpha1_beats"] if beats is None else beats
    x = np.asarray(x, dtype=float)
    _require_variability(x, resolution_ms)
    y = np.cumsum(x - x.mean())
    ns, fn = [], []
    for n in range(lo, hi + 1):
        k = len(y) // n
        if k < 2:
            break
        boxes = y[:k * n].reshape(k, n)
        t = np.arange(n)
        slope, icpt = np.polyfit(t, boxes.T, 1)              # one line per box
        resid = boxes - (np.outer(slope, t) + icpt[:, None])
        ns.append(n)
        fn.append(np.sqrt(np.mean(resid ** 2)))
    if len(ns) < 3 or min(fn) <= 0:
        raise NotApplicable("series too short for DFA alpha1")
    return float(np.polyfit(np.log(ns), np.log(fn), 1)[0])


def _rr(ctx) -> np.ndarray:
    rr = rr_ms(ctx.r_peaks, ctx.fs)
    if len(rr) < 2:
        raise NotApplicable("fewer than 3 detected beats")
    return rr


def _res(ctx) -> float:
    return 1000.0 / ctx.fs


# ------------------------------------------------------ on every record ---
@feature(name="n_beats", family="rate", units="beats", leads=["any"], requires=["r_peaks"],
         desc="Number of detected QRS complexes")
def n_beats(ctx):
    return int(len(ctx.r_peaks))


@feature(name="rr_intervals_ms", family="rate", units="ms", leads=["any"],
         requires=["r_peaks"], min_duration_s=10,
         desc="RR interval series between consecutive detected beats")
def rr_intervals_ms(ctx):
    return _rr(ctx).round(1).tolist()


@feature(name="mean_hr_bpm", family="rate", units="bpm", leads=["any"],
         requires=["r_peaks"], min_duration_s=10,
         desc="Mean heart rate, 60000 / mean RR interval")
def mean_hr_bpm(ctx):
    return float(60000.0 / np.mean(_rr(ctx)))


@feature(name="rmssd_ms", family="rate", units="ms", leads=["any"], requires=["r_peaks"],
         min_duration_s=10, cite="taskforce1996hrv",
         desc="Root mean square of successive RR differences")
def rmssd_ms(ctx):
    return rmssd(_rr(ctx))


# ------------------------------------------- long recordings only (Rule 5) ---
@feature(name="sdnn_ms", family="rate", units="ms", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="taskforce1996hrv",
         desc="Standard deviation of RR intervals (5-min definition)")
def sdnn_ms(ctx):
    return sdnn(_rr(ctx))


@feature(name="pnn50_pct", family="rate", units="%", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="taskforce1996hrv",
         desc="Share of successive RR differences > 50 ms")
def pnn50_pct(ctx):
    return pnn50(_rr(ctx))


@feature(name="lf_ms2", family="rate", units="ms²", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="taskforce1996hrv",
         desc="Low-frequency RR power, 0.04-0.15 Hz")
def lf_ms2(ctx):
    return lf_hf(_rr(ctx), _res(ctx))[0]


@feature(name="hf_ms2", family="rate", units="ms²", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="taskforce1996hrv",
         desc="High-frequency RR power, 0.15-0.40 Hz")
def hf_ms2(ctx):
    return lf_hf(_rr(ctx), _res(ctx))[1]


@feature(name="lf_hf_ratio", family="rate", units="", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="taskforce1996hrv", desc="LF / HF power ratio")
def lf_hf_ratio(ctx):
    return lf_hf(_rr(ctx), _res(ctx))[2]


@feature(name="sampen", family="rate", units="", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="richman2000sampen",
         desc="Sample entropy of the RR series (m = 2, r = 0.2 SD)")
def sampen_feature(ctx):
    return sampen(_rr(ctx), resolution_ms=_res(ctx))


@feature(name="dfa_alpha1", family="rate", units="", leads=["any"], requires=["r_peaks"],
         min_duration_s=LONG_S, cite="peng1995dfa",
         desc="Short-term detrended fluctuation exponent, 4-16 beats")
def dfa_alpha1_feature(ctx):
    return dfa_alpha1(_rr(ctx), resolution_ms=_res(ctx))
