"""
Wavelet-based single-lead delineator of Martínez et al. (IEEE TBME 2004;
51(4):570-581, cite key martinez2004wavedet), plus the multilead selection rule
they applied to the 12-lead CSE database (Laguna et al. 1994).

Thresholds are the paper's Appendix (A.3-A.14) and live in config.yaml under
`wavedet:`. Where the paper leaves a choice open (search windows) the value is
in config.yaml too and marked as ours.

Deviations from the paper, stated for the Methods section:
  * QRS detection: our multilead R-peak detector (delineate.detect_r_peaks) is
    used instead of the paper's modulus-maximum-line detector; every lead is
    delineated around the same R positions.
  * Sampling rate: the paper resamples its 250 Hz filter bank to other rates.
    We keep the same quadratic-spline filters and shift the scale index by
    log2(fs / 250), which gives the same analogue bands at 500 Hz exactly
    (scale 2^k at 250 Hz == scale 2^(k+1) at 500 Hz).
  * A lead's P wave is kept only if its signal deflection is >= bounds.p_min_mv,
    and a global mark needs k + 1 leads to report the wave (ours).
  * Significant slopes of a P or T wave must lie within p_span_ms / t_span_ms of
    the dominant slope (ours; the paper relies on wave-sized search windows,
    which it does not specify).
  * Significance thresholds (gamma) are relative to the largest local modulus
    maximum inside the window, not the window-wide maximum, so the monotonic
    leakage of an adjacent QRS into a P or T window cannot set the threshold.
Scale numbers below are the paper's (250 Hz) numbering.
"""

import math

import numpy as np

from . import config

# Quadratic spline (Mallat & Zhong; Li et al. 1995): h = [1,3,3,1]/8, g = 2[1,-1]
_H = np.array([1, 3, 3, 1]) / 8.0
_G = np.array([1, -1]) * 2.0


def _upsample(f, j):
    """Insert 2^j - 1 zeros between taps (algorithme a trous)."""
    if j == 0:
        return f
    out = np.zeros((len(f) - 1) * 2 ** j + 1)
    out[:: 2 ** j] = f
    return out


def dyadic_wt(x: np.ndarray, fs: int, n_scales: int = 5) -> dict[int, np.ndarray]:
    """{k: W_{2^k} x} for k = 1..n_scales in the paper's 250 Hz numbering,
    each aligned to x (zero-phase up to half a sample)."""
    shift = int(round(math.log2(fs / 250)))
    if shift < 0 or abs(fs / 250 - 2 ** shift) > 1e-9:
        raise ValueError(f"wavedet supports fs = 250 * 2^n Hz, got {fs}")
    x = np.asarray(x, dtype=float)
    out, low = {}, np.array([1.0])
    for j in range(n_scales + shift):
        psi = np.convolve(low, _upsample(_G, j))           # equivalent wavelet filter
        k = j + 1 - shift
        if k >= 1:
            delay = (len(psi) - 1) / 2
            y = np.convolve(x, psi, mode="full")
            d = int(math.floor(delay))
            out[k] = y[d:d + len(x)]
        low = np.convolve(low, _upsample(_H, j))
    return out


# ------------------------------------------------------------ helpers -------
def _local_maxima(a, lo, hi):
    """Indices of local maxima of |a| strictly inside [lo, hi)."""
    lo, hi = max(lo, 1), min(hi, len(a) - 1)
    if hi - lo < 1:
        return np.array([], dtype=int)
    m = np.abs(a[lo - 1:hi + 1])
    i = np.where((m[1:-1] > m[:-2]) & (m[1:-1] >= m[2:]))[0]
    return i + lo


def _boundary(w, n, xi, step, limit):
    """From n, walk by `step` until |w| < xi or |w| has a local minimum.
    The paper takes whichever comes first (nearest the wave)."""
    a = np.abs(w)
    i = n
    while (i - n) * step < limit:
        i += step
        if not 0 < i < len(a) - 1:
            return math.nan
        if a[i] < xi:
            return float(i)
        if a[i] <= a[i - 1] and a[i] <= a[i + 1]:
            return float(i)
    return math.nan


def _near_dominant(w, sig, span):
    """Keep significant maxima within `span` samples of the dominant one, so a
    wiggle elsewhere in the window is not taken as the wave's first/last slope."""
    if len(sig) == 0:
        return sig
    top = sig[np.argmax(np.abs(w[sig]))]
    return sig[np.abs(sig - top) <= span]


def _zero_cross(w, a, b):
    lo, hi = sorted((a, b))
    s = np.sign(w[lo:hi + 1])
    z = np.where(s[:-1] * s[1:] <= 0)[0]
    return float(lo + z[0]) if len(z) else math.nan


# ------------------------------------------------------------ one lead ------
def delineate_lead(x, fs, r_peaks, cfg=None) -> dict[str, np.ndarray]:
    """Per-beat fiducials (float samples, NaN if absent) for one lead."""
    cfg = cfg or config.load()
    c = cfg["wavedet"]
    p_min_mv = cfg["bounds"]["p_min_mv"]
    x = np.asarray(x, dtype=float)
    ms = lambda v: int(round(v * fs / 1000))                     # noqa: E731
    W = dyadic_wt(x, fs)
    w2, w4, w5 = W[2], W[4], W[5]
    n = len(r_peaks)
    keys = ("p_onset", "p_peak", "p_offset", "qrs_onset", "qrs_offset", "t_peak", "t_offset")
    out = {k: np.full(n, np.nan) for k in keys}
    rr = np.diff(r_peaks)
    rr_med = float(np.median(rr)) if len(rr) else fs

    for b, r in enumerate(r_peaks):
        rr_b = float(rr[b - 1]) if b > 0 else rr_med
        rr_next = float(rr[b]) if b < len(rr) else rr_med
        # --- QRS (scale 2^2) -----------------------------------------------
        lo, hi = r + ms(c["qrs_window_ms"][0]), r + ms(c["qrs_window_ms"][1])
        mx = _local_maxima(w2, lo, hi)
        if len(mx) == 0:
            continue
        top = np.abs(w2[mx]).max()
        pre = mx[(mx <= r) & (np.abs(w2[mx]) > c["gamma_qrs_pre"] * top)]
        post = mx[(mx >= r) & (np.abs(w2[mx]) > c["gamma_qrs_post"] * top)]
        sig = np.concatenate([pre, post])
        if len(sig) == 0:
            continue
        n_first, n_last = int(sig.min()), int(sig.max())
        xi_on = (c["xi_qrs_on"][0] if w2[n_first] > 0 else c["xi_qrs_on"][1]) * abs(w2[n_first])
        xi_end = (c["xi_qrs_end"][0] if w2[n_last] > 0 else c["xi_qrs_end"][1]) * abs(w2[n_last])
        q_on = _boundary(w2, n_first, xi_on, -1, ms(c["boundary_search_ms"]))
        q_off = _boundary(w2, n_last, xi_end, +1, ms(c["boundary_search_ms"]))
        out["qrs_onset"][b], out["qrs_offset"][b] = q_on, q_off

        seg_rms = lambda w, a, z: float(np.sqrt(np.mean(w[max(a, 0):min(z, len(w))] ** 2)))  # noqa: E731

        # --- T wave (scale 2^4, else 2^5) -----------------------------------
        t_lo = int(q_off if not math.isnan(q_off) else r) + ms(c["t_window_ms"][0])
        t_hi = r + min(int(c["t_window_rr_frac"] * rr_next), ms(c["t_window_ms"][1]))
        for w in (w4, w5):
            mx = _local_maxima(w, t_lo, t_hi)
            eps = c["eps_t"] * seg_rms(w, r, r + int(rr_next))
            if (np.abs(w[mx]) > eps).sum() < 2:
                continue
            top = np.abs(w[mx]).max()
            sig = _near_dominant(w, mx[np.abs(w[mx]) > c["gamma_t"] * top], ms(c["t_span_ms"]))
            if len(sig) < 1:
                continue
            nf, nl = int(sig.min()), int(sig.max())
            big = sig[np.argsort(-np.abs(w[sig]))[:2]]
            if len(big) == 2 and np.sign(w[big[0]]) != np.sign(w[big[1]]):
                out["t_peak"][b] = _zero_cross(w, big[0], big[1])
            out["t_offset"][b] = _boundary(w, nl, c["xi_t_end"] * abs(w[nl]), +1,
                                           ms(c["boundary_search_ms"]) * 2)
            break

        # --- P wave (scale 2^4) ----------------------------------------------
        if math.isnan(q_on):
            continue
        p_hi = int(q_on)
        p_lo = r - min(int(c["p_window_rr_frac"] * rr_b), ms(c["p_window_ms"]))
        if b > 0 and not math.isnan(out["t_offset"][b - 1]):
            p_lo = max(p_lo, int(out["t_offset"][b - 1]))
        mx = _local_maxima(w4, p_lo, p_hi)
        eps = c["eps_p"] * seg_rms(w4, r - int(rr_b), r)
        if (np.abs(w4[mx]) > eps).sum() < 2:
            continue
        top = np.abs(w4[mx]).max()
        sig = _near_dominant(w4, mx[np.abs(w4[mx]) > c["gamma_p"] * top], ms(c["p_span_ms"]))
        if len(sig) < 2:
            continue
        nf, nl = int(sig.min()), int(sig.max())
        big = sig[np.argsort(-np.abs(w4[sig]))[:2]]
        if np.sign(w4[big[0]]) != np.sign(w4[big[1]]):
            out["p_peak"][b] = _zero_cross(w4, big[0], big[1])
        p_on = _boundary(w4, nf, c["xi_p_on"] * abs(w4[nf]), -1, ms(c["boundary_search_ms"]))
        p_off = _boundary(w4, nl, c["xi_p_end"] * abs(w4[nl]), +1, ms(c["boundary_search_ms"]))
        if math.isnan(p_on) or math.isnan(p_off) or p_off <= p_on:
            continue
        # eps_P (A.4) is permissive enough to accept QRS/T leakage when there is
        # no P wave at all; require a P deflection of at least p_min_mv (ours).
        seg = x[int(p_on):int(p_off) + 1]
        if np.max(np.abs(seg - x[int(p_on)])) < p_min_mv:
            out["p_peak"][b] = math.nan
            continue
        out["p_onset"][b], out["p_offset"][b] = p_on, p_off
    return out


# ------------------------------------------------------------ multilead -----
def select_multilead(marks: np.ndarray, onset: bool, k: int, delta: float) -> float:
    """Laguna et al. 1994 rule: sort the single-lead marks; the onset (end) is the
    first (last) mark whose k nearest neighbours lie within `delta` samples.
    Falls back to the median when no mark qualifies."""
    m = np.sort(marks[~np.isnan(marks)])
    if len(m) == 0:
        return math.nan
    if not onset:
        m = m[::-1]
    for i in range(len(m)):
        d = np.sort(np.abs(m - m[i]))[1:k + 1]
        if len(d) == k and d.max() <= delta:
            return float(m[i])
    return float(np.median(m))


def delineate_multilead(sig: np.ndarray, fs: int, r_peaks, cfg=None, return_per_lead=False):
    """sig: (n_samples, n_leads). Returns global per-beat fiducials, and with
    return_per_lead=True also the single-lead marks they were pooled from
    (list, one dict per column of sig). The pooled marks are the same either way."""
    cfg = cfg or config.load()
    c = cfg["wavedet"]
    per = [delineate_lead(sig[:, j], fs, r_peaks, cfg) for j in range(sig.shape[1])]
    out = {}
    for key in per[0]:
        a = np.vstack([p[key] for p in per])                   # (leads, beats)
        if key.endswith("_peak"):
            out[key] = np.nanmedian(a, axis=0) if a.size else a
            continue
        rule = c["multilead"]
        onset = key.endswith("onset")
        # a global mark needs k + 1 leads agreeing it exists (ours), so one lead's
        # spurious wave cannot become the recording's wave
        k = rule["k"][key] if isinstance(rule["k"], dict) else rule["k"]
        if sig.shape[1] > k:
            too_few = np.sum(~np.isnan(a), axis=0) < k + 1
            a = np.where(too_few[None, :], np.nan, a)
        if sig.shape[1] == 1 or rule["method"] == "none":
            out[key] = a[0]
        elif rule["method"] == "extreme":
            out[key] = np.nanmin(a, 0) if onset else np.nanmax(a, 0)
        else:
            delta = rule["delta_ms"][key] * fs / 1000
            out[key] = np.array([select_multilead(a[:, i], onset, k, delta)
                                 for i in range(a.shape[1])])
    return (out, per) if return_per_lead else out
