"""Notch + band-pass, zero phase. Shared by every source, so PTB-XL and the
device are filtered identically."""

from functools import lru_cache

import numpy as np
from scipy import signal as ss

from . import config
from .recording import Recording


@lru_cache(maxsize=None)
def _design(fs: int, notch_hz: float, notch_q: float, band: tuple, order: int):
    lo, hi = band
    nyq = fs / 2
    if not 0 < lo < hi < nyq:
        raise ValueError(f"band {band} Hz is not valid at fs={fs} Hz")
    sos = ss.butter(order, [lo, hi], btype="bandpass", fs=fs, output="sos")
    notch = None
    if notch_hz < nyq:
        b, a = ss.iirnotch(notch_hz, notch_q, fs=fs)
        notch = ss.tf2sos(b, a)
    return sos, notch


def filter_signal(x: np.ndarray, fs: int, cfg: dict | None = None) -> np.ndarray:
    """x: (n_samples, n_leads) in mV -> filtered copy, same shape."""
    f = (cfg or config.load())["filter"]
    sos, notch = _design(int(fs), float(f["notch_hz"]), float(f["notch_q"]),
                         tuple(f["band"]), int(f["order"]))
    y = np.asarray(x, dtype=float) - np.median(x, axis=0)   # remove DC before filtfilt
    if notch is not None:
        y = ss.sosfiltfilt(notch, y, axis=0)
    return ss.sosfiltfilt(sos, y, axis=0)


def preprocess(rec: Recording, cfg: dict | None = None) -> Recording:
    return Recording(filter_signal(rec.signal, rec.fs, cfg), rec.fs, rec.leads,
                     rec.duration_s, rec.source, rec.record_id,
                     {**rec.meta, "filtered": True})
