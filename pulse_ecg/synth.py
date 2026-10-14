"""
Synthetic 12-lead ECG with KNOWN answers — the ground truth for tests/.

Every wave is a half-sine (or a chain of half-sines) with exact support, so
its onset and offset are known to the sample. Limb leads are projections of a
frontal-plane axis onto the hexaxial reference, so they obey Einthoven's law
(III = II - I) and the P axis has a known true value. It is not physiological
in detail; it only has to make "did we recover the value we put in" testable.
"""

import numpy as np

from .recording import STANDARD_LEADS, Recording

# Precordial gains per wave: (P, QRS, T). V1 QRS is an rS complex, so negative.
_PRECORDIAL = {"V1": (0.5, -0.6, 0.2), "V2": (0.6, 0.4, 0.8), "V3": (0.6, 0.8, 0.9),
               "V4": (0.5, 1.2, 0.9), "V5": (0.4, 1.1, 0.7), "V6": (0.4, 0.9, 0.6)}


def _add_half_sine(col, t, onset, width, amp):
    """col[support] += amp x half-sine, in place. t is uniform and sorted, so the
    support is found by bisection rather than a full-length mask per wave."""
    lo, hi = np.searchsorted(t, onset, "left"), np.searchsorted(t, onset + width, "left")
    col[lo:hi] += amp * np.sin(np.pi * (t[lo:hi] - onset) / width)


def _gain(lead, wave_idx, axis_deg):
    """Lead gain for a wave whose frontal axis is axis_deg.

    Only I and II are projections (at 0 and 60 degrees); the other four limb
    leads are derived from them exactly as an electrocardiograph derives them
    (Einthoven and Goldberger). Projecting aVR/aVL/aVF directly would make them
    sqrt(3)/2 too small. The frontal axis is then recovered as
    atan2(2*aVF/sqrt(3), I)."""
    th = np.deg2rad(axis_deg)
    i, ii = np.cos(th), np.cos(th - np.pi / 3)
    limb = {"I": i, "II": ii, "III": ii - i, "aVR": -(i + ii) / 2,
            "aVL": i - ii / 2, "aVF": ii - i / 2}
    if lead in limb:
        return float(limb[lead])
    return _PRECORDIAL[lead][wave_idx]


def synth(*, hr=60.0, duration_s=10.0, fs=500, leads=STANDARD_LEADS,
          p_duration_ms=100.0, pr_ms=160.0, qrs_duration_ms=90.0, qt_ms=400.0,
          p_amp_mv=0.15, r_amp_mv=1.0, t_amp_mv=0.3,
          p_axis_deg=50.0, qrs_axis_deg=45.0, t_axis_deg=45.0,
          biphasic_v1=False, v1_terminal_frac=0.5, v1_terminal_amp_mv=0.1,
          noise_mv=0.0, seed=0, return_truth=False):
    """
    Build a synthetic Recording. PR is P onset to QRS onset; QT is QRS onset to
    T offset. With biphasic_v1=True, lead V1's P wave is positive for its first
    (1 - v1_terminal_frac) and negative, with amplitude v1_terminal_amp_mv, for
    the rest; its PTFV1 is then -(amplitude in µV) x (duration in ms).
    """
    if pr_ms <= p_duration_ms:
        raise ValueError("pr_ms must exceed p_duration_ms (P ends before QRS starts)")
    if qt_ms <= qrs_duration_ms:
        raise ValueError("qt_ms must exceed qrs_duration_ms")
    if not 0 < v1_terminal_frac < 1:
        raise ValueError("v1_terminal_frac must be in (0, 1)")

    all_leads = STANDARD_LEADS
    n = int(round(duration_s * fs))
    t = np.arange(n) / fs
    rr = 60.0 / hr
    p_w, pr, qrs_w, qt = (x / 1000 for x in (p_duration_ms, pr_ms, qrs_duration_ms, qt_ms))
    t_w = 0.5 * (qt - qrs_w)

    first = pr + 0.2
    onsets = []
    k = 0
    while first + k * rr + qt < duration_s - 0.1:
        onsets.append(first + k * rr)
        k += 1
    if not onsets:
        raise ValueError("recording too short for one complete beat")

    sig = np.zeros((n, len(all_leads)))
    for j, lead in enumerate(all_leads):
        gp, gq, gt = (_gain(lead, 0, p_axis_deg), _gain(lead, 1, qrs_axis_deg),
                      _gain(lead, 2, t_axis_deg))
        for q0 in onsets:
            p0 = q0 - pr
            if lead == "V1" and biphasic_v1:
                w1 = p_w * (1 - v1_terminal_frac)
                _add_half_sine(sig[:, j], t, p0, w1, p_amp_mv * gp)
                _add_half_sine(sig[:, j], t, p0 + w1, p_w - w1, -v1_terminal_amp_mv)
            else:
                _add_half_sine(sig[:, j], t, p0, p_w, p_amp_mv * gp)
            # QRS: Q (15 %), R (50 %), S (35 %) of the width, contiguous
            _add_half_sine(sig[:, j], t, q0, 0.15 * qrs_w, -0.10 * r_amp_mv * gq)
            _add_half_sine(sig[:, j], t, q0 + 0.15 * qrs_w, 0.50 * qrs_w, r_amp_mv * gq)
            _add_half_sine(sig[:, j], t, q0 + 0.65 * qrs_w, 0.35 * qrs_w, -0.25 * r_amp_mv * gq)
            _add_half_sine(sig[:, j], t, q0 + qt - t_w, t_w, t_amp_mv * gt)

    if noise_mv > 0:
        sig += np.random.default_rng(seed).normal(0.0, noise_mv, sig.shape)

    rec = Recording.from_array(sig, fs, all_leads, "synthetic", f"synth-{seed}")
    if tuple(leads) != all_leads:
        rec = rec.select(leads)
    if not return_truth:
        return rec

    q = np.asarray(onsets)
    idx = lambda secs: np.round(secs * fs).astype(int)  # noqa: E731
    truth = {
        "p_onset": idx(q - pr), "p_offset": idx(q - pr + p_w),
        "qrs_onset": idx(q), "qrs_offset": idx(q + qrs_w), "t_offset": idx(q + qt),
        "n_beats": len(q), "mean_hr_bpm": hr,
        "p_duration_ms": p_duration_ms, "pr_ms": pr_ms,
        "qrs_duration_ms": qrs_duration_ms, "qt_ms": qt_ms,
        "p_axis_deg": p_axis_deg,
        "ptfv1_uv_ms": (-v1_terminal_amp_mv * 1000 * p_duration_ms * v1_terminal_frac
                        if biphasic_v1 else 0.0),
    }
    return rec, truth
