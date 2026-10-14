"""
Delineation: R-peaks -> median beat -> vector magnitude -> wave boundaries.

Rule 3: global intervals are measured on the vector magnitude of every lead
present (12 for PTB-XL, 1 for the device), never on a single chosen lead.

    d = delineate(preprocess(rec))
    d.points        median-beat indices (float samples): p_onset, p_peak, p_offset,
                    qrs_onset, r_peak, qrs_offset, t_peak, t_offset  (NaN if absent)
    d.r_peaks       R-peak sample index of every detected beat in the recording
    d.beat_points() the median-beat points mapped onto every beat (for Gate 2)

Methods (every constant in config.yaml):
  QRS   spatial velocity (norm of the per-lead derivative) falls to within
        bounds.qrs_frac of its isoelectric floor (floor + frac x (peak - floor))
        for bounds.min_run_ms. The floor is the bounds.floor_pct percentile of the
        median beat's spatial velocity.
  P, T  vector magnitude falls to its local floor + bounds.p_frac / t_frac x the
        wave's height and stays there for min_run_ms. The floor is the local
        minimum, because the VM of 12 noisy leads never reaches zero.
"""

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import signal as ss

from . import config
from .recording import Recording

POINTS = ("p_onset", "p_peak", "p_offset", "qrs_onset", "r_peak",
          "qrs_offset", "t_peak", "t_offset")


@dataclass
class Delineation:
    fs: int
    r_peaks: np.ndarray                    # absolute samples, one per detected beat
    r0: int                                # R position inside the median beat
    median_beat: np.ndarray | None         # (win, n_leads), baseline-corrected
    vm: np.ndarray | None                  # vector magnitude of median_beat
    points: dict = field(default_factory=dict)
    beat_template_corr: float = math.nan
    n_beats_median: int = 0                # beats complete enough to enter the median
    baseline: np.ndarray | None = None     # per-lead value subtracted from median_beat
    per_beat: dict | None = None           # native per-beat fiducials (wavedet)
    per_lead: dict | None = None           # {lead: {point: median-beat sample}} single-lead marks
    method: str = "threshold"
    reason: str = ""                       # why delineation failed, if it did

    @property
    def ok(self) -> bool:
        return not self.reason

    def beat_points(self) -> dict[str, np.ndarray]:
        """Fiducials of every beat: native per-beat marks when the method gives
        them (wavedet), else the median-beat points placed on each R-peak."""
        if self.per_beat is not None:
            return {**self.per_beat, "r_peak": self.r_peaks.astype(float)}
        return {k: self.r_peaks + (v - self.r0) for k, v in self.points.items()}

    def context_extra(self) -> dict:
        return {"median_beat": self.median_beat, "vm": self.vm, "r_peaks": self.r_peaks,
                "r0": self.r0, "beat_template_corr": self.beat_template_corr,
                "per_beat": self.per_beat, "per_lead": self.per_lead}


# ------------------------------------------------------------------ R-peaks ---
def detect_r_peaks(x: np.ndarray, fs: int, cfg: dict | None = None) -> np.ndarray:
    """Energy of the multi-lead derivative, integrated, then peak-picked.
    x: filtered (n_samples, n_leads)."""
    c = (cfg or config.load())["rpeak"]
    d = np.diff(x, axis=0, prepend=x[:1])
    e = np.sum(d ** 2, axis=1)
    w = max(1, int(round(c["integ_ms"] * fs / 1000)))
    e = np.convolve(e, np.ones(w) / w, mode="same")
    dist = max(1, int(round(c["refractory_ms"] * fs / 1000)))
    cand, props = ss.find_peaks(e, distance=dist, height=0)
    if len(cand) == 0:
        return cand
    ref = np.percentile(props["peak_heights"], 90)
    cand = cand[props["peak_heights"] >= c["thresh_frac"] * ref]

    # refine to the largest absolute deflection nearby
    mag = np.sqrt(np.sum(x ** 2, axis=1))
    half = int(round(c["refine_ms"] * fs / 1000))
    out = []
    for p in cand:
        lo, hi = max(0, p - half), min(len(mag), p + half + 1)
        out.append(lo + int(np.argmax(mag[lo:hi])))
    out = np.unique(out)
    keep = [out[0]] if len(out) else []           # re-apply refractory after refining
    for p in out[1:]:
        if p - keep[-1] >= dist:
            keep.append(p)
        elif mag[p] > mag[keep[-1]]:
            keep[-1] = p
    return np.asarray(keep, dtype=int)


# ------------------------------------------------------------- median beat ---
def median_beat(x: np.ndarray, r_peaks: np.ndarray, fs: int, cfg: dict | None = None):
    """-> (median beat (win, n_leads), r0, beats used, median beat-template correlation)."""
    b = (cfg or config.load())["beat"]
    pre, post = int(round(b["pre_s"] * fs)), int(round(b["post_s"] * fs))
    used = np.asarray([r for r in r_peaks if r - pre >= 0 and r + post <= len(x)], dtype=int)
    if len(used) < b["min_beats"]:
        return None, pre, used, math.nan
    beats = np.stack([x[r - pre:r + post] for r in used])       # (k, win, leads)
    med = np.median(beats, axis=0)
    flat = med.ravel() - med.mean()
    corr = [np.corrcoef(flat, bt.ravel() - bt.mean())[0, 1] for bt in beats]
    return med, pre, used, float(np.nanmedian(corr))


# --------------------------------------------------------------- helpers ---
def _quiet_edge(a, start, thr, run, step, stop):
    """Walk from `start` in direction `step` (+1/-1) to the first index i such that
    a[i], a[i+step], ... (run samples) are all < thr. Returns i or None."""
    i = start
    while (step < 0 and i - run + 1 >= stop) or (step > 0 and i + run - 1 <= stop):
        seg = a[i - run + 1:i + 1] if step < 0 else a[i:i + run]
        if np.all(seg < thr):
            return i
        i += step
    return None


def _wave_edge(vm, peak, frac, run, step, stop):
    """Edge of the wave peaking at `peak`: the first sample, walking away from the
    peak, where vm falls to floor + frac x (vm[peak] - floor) and stays there for
    `run` samples. The floor is the local minimum of vm between peak and stop.
    Returns a float sample index, or NaN."""
    seg = vm[stop:peak + 1] if step < 0 else vm[peak:stop + 1]
    if len(seg) < run + 1:
        return math.nan
    floor = seg.min()
    i = _quiet_edge(vm, peak, floor + frac * (vm[peak] - floor), run, step, stop)
    return math.nan if i is None else float(i)


def _median_rel(marks, r_peaks) -> float:
    """Median position of per-beat marks relative to their R-peak (NaN if none)."""
    return float(np.nanmedian(marks - r_peaks)) if np.any(~np.isnan(marks)) else math.nan


# -------------------------------------------------------------- delineate ---
def delineate(rec: Recording, cfg: dict | None = None) -> Delineation:
    """`rec` should already be preprocessed (preprocess.preprocess).
    config `delineation.method`: "wavedet" (Martínez et al. 2004, per lead, then
    the multilead rule; the validated default) or "threshold" (median-beat vector
    magnitude; failed Gate 2, kept for comparison)."""
    cfg = cfg or config.load()
    d = _delineate_threshold(rec, cfg)
    if cfg["delineation"]["method"] == "threshold" or d.median_beat is None:
        return d
    if cfg["delineation"]["method"] != "wavedet":
        raise ValueError(f"unknown delineation.method {cfg['delineation']['method']!r}")
    import warnings
    from .wavedet import delineate_multilead
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # all-NaN beats
        per, singles = delineate_multilead(rec.signal, rec.fs, d.r_peaks, cfg,
                                           return_per_lead=True)
        rel = {k: _median_rel(v, d.r_peaks) for k, v in per.items()}
        d.per_lead = {lead: {k: d.r0 + _median_rel(v, d.r_peaks) for k, v in m.items()}
                      for lead, m in zip(rec.leads, singles)}
    d.per_beat, d.method = per, "wavedet"
    d.points = {k: math.nan for k in POINTS}
    d.points.update({k: d.r0 + v for k, v in rel.items() if k in POINTS}, r_peak=float(d.r0))
    if d.median_beat is not None and not math.isnan(d.points["qrs_onset"]):
        # re-reference the median beat to the isoelectric run before wavedet's QRS onset
        run = max(2, int(round(cfg["bounds"]["min_run_ms"] * rec.fs / 1000)))
        q = int(round(d.points["qrs_onset"]))
        raw = d.median_beat + d.baseline
        d.baseline = raw[max(0, q - run):max(q, 1)].mean(axis=0)
        d.median_beat = raw - d.baseline
        d.vm = np.sqrt(np.sum(d.median_beat ** 2, axis=1))
    if all(math.isnan(d.points[k]) for k in ("qrs_onset", "qrs_offset")):
        d.reason = "wavedet found no QRS boundaries"
    else:
        d.reason = ""
    return d


def _delineate_threshold(rec: Recording, cfg: dict) -> Delineation:
    fs, x = rec.fs, rec.signal
    ms = lambda v: int(round(v * fs / 1000))                     # noqa: E731
    bd = cfg["bounds"]
    run = max(2, ms(bd["min_run_ms"]))
    nan_points = {k: math.nan for k in POINTS}

    r_peaks = detect_r_peaks(x, fs, cfg)
    med, r0, used, corr = median_beat(x, r_peaks, fs, cfg)
    d = Delineation(fs, r_peaks, r0, None, None, dict(nan_points), corr, len(used))
    if med is None:
        d.reason = f"fewer than {cfg['beat']['min_beats']} complete beats"
        return d
    d.median_beat, d.baseline = med, np.zeros(med.shape[1])   # replaced below if QRS found

    # --- QRS on spatial velocity (baseline-independent) --------------------
    sv = np.sqrt(np.sum(np.gradient(med, axis=0) ** 2, axis=1))
    lo, hi = max(0, r0 - ms(60)), min(len(sv), r0 + ms(60) + 1)
    floor = np.percentile(sv, bd["floor_pct"])        # isoelectric level of the beat
    thr = floor + bd["qrs_frac"] * (sv[lo:hi].max() - floor)
    on = _quiet_edge(sv, r0, thr, run, -1, 0)
    off = _quiet_edge(sv, r0, thr, run, +1, len(sv) - 1)
    if on is None or off is None:
        d.reason = "QRS boundaries not found"
        return d
    qrs_on, qrs_off = on + 1, off - 1

    # --- baseline: the isoelectric run just before QRS onset ---------------
    base = med[max(0, qrs_on - run):qrs_on].mean(axis=0)
    mc = med - base
    vm = np.sqrt(np.sum(mc ** 2, axis=1))
    d.median_beat, d.vm, d.baseline = mc, vm, base
    pts = d.points
    pts.update(qrs_onset=float(qrs_on), qrs_offset=float(qrs_off),
               r_peak=float(lo + np.argmax(vm[lo:hi])))

    # --- T wave: peak, then return to the local floor -----------------------
    t_lo = qrs_off + ms(bd["t_search_ms"][0])
    t_hi = min(len(vm) - 1, qrs_on + ms(bd["t_search_ms"][1]))
    if t_hi - t_lo > run:
        tp = t_lo + int(np.argmax(vm[t_lo:t_hi]))
        t_off = _wave_edge(vm, tp, bd["t_frac"], run, +1, min(len(vm) - 1, tp + ms(250)))
        if not math.isnan(t_off):
            pts.update(t_peak=float(tp), t_offset=t_off)

    # --- P wave: before QRS onset, after the previous beat's T ------------
    p_hi = qrs_on - run
    p_lo = max(0, qrs_on - ms(bd["p_search_ms"]))
    if len(d.r_peaks) > 1 and not math.isnan(pts["t_offset"]):
        rr = int(np.median(np.diff(d.r_peaks)))   # previous T ends ~rr before this one
        p_lo = max(p_lo, int(math.ceil(pts["t_offset"])) - rr + run)
    if p_hi - p_lo > run:
        pp = p_lo + int(np.argmax(vm[p_lo:p_hi]))
        amp = vm[pp] - vm[p_lo:p_hi].min()
        if amp >= bd["p_min_mv"] and pp not in (p_lo, p_hi - 1):
            p_on = _wave_edge(vm, pp, bd["p_frac"], run, -1, max(0, pp - ms(bd["p_search_ms"])))
            p_off = _wave_edge(vm, pp, bd["p_frac"], run, +1, qrs_on)
            if not (math.isnan(p_on) or math.isnan(p_off)):
                pts.update(p_onset=p_on, p_peak=float(pp), p_offset=p_off)
    return d
