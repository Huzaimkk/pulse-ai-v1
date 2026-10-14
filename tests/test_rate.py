"""Rate family: known answers on synthetic R-peak trains and RR series.

Long-window HRV (SDNN, pNN50, LF/HF, SampEn, DFA alpha1) must stay dormant on a
10 s record (Rule 5) and fire on a device-length one — same code, gated by
min_duration_s in the declaration."""

import math

import numpy as np
import pytest

from pulse_ecg.extract import extract
from pulse_ecg.features import rate
from pulse_ecg.synth import synth

LONG = ("sdnn_ms", "pnn50_pct", "lf_ms2", "hf_ms2", "lf_hf_ratio", "sampen", "dfa_alpha1")


@pytest.mark.parametrize("hr", [45, 60, 90])
def test_mean_hr_and_rr_series_recover_synthetic_rate(hr):
    out = extract(synth(hr=hr), strict=True)
    assert out["mean_hr_bpm"] == pytest.approx(hr, abs=1.0)
    rr = np.asarray(out["rr_intervals_ms"])
    assert len(rr) == out["n_beats"] - 1
    assert np.all(np.abs(rr - 60000 / hr) <= 4)          # 2 samples at 500 Hz
    assert out["rmssd_ms"] < 5                            # perfectly regular rhythm


def test_rmssd_sdnn_pnn50_match_their_definitions():
    rr = np.array([800.0, 810, 790, 800])
    assert rate.rmssd(rr) == pytest.approx(math.sqrt((100 + 400 + 100) / 3))
    assert rate.sdnn(rr) == pytest.approx(np.std(rr, ddof=1))
    assert rate.pnn50(np.array([800.0, 860, 800, 830])) == pytest.approx(200 / 3)


def _tachogram(f_hz, seconds=300, mean_ms=1000.0, amp_ms=50.0):
    t, rr = 0.0, []
    while t < seconds:
        rr.append(mean_ms + amp_ms * math.sin(2 * math.pi * f_hz * t))
        t += rr[-1] / 1000
    return np.asarray(rr)


def test_lf_hf_ratio_follows_the_modulation_frequency():
    lf, hf, ratio = rate.lf_hf(_tachogram(0.10))          # Task Force LF band
    assert ratio > 5 and lf > hf
    lf, hf, ratio = rate.lf_hf(_tachogram(0.25))          # HF band
    assert ratio < 0.2


def test_sample_entropy_is_high_for_noise_and_low_for_a_regular_series():
    rng = np.random.default_rng(0)
    assert 1.9 < rate.sampen(rng.normal(size=1000)) < 2.5   # Gaussian white noise ~2.2
    assert rate.sampen(_tachogram(0.1)) < 0.5


def test_dfa_alpha1_separates_white_noise_from_a_random_walk():
    rng = np.random.default_rng(1)
    noise = rng.normal(size=2000)
    assert rate.dfa_alpha1(noise) == pytest.approx(0.5, abs=0.15)
    assert rate.dfa_alpha1(np.cumsum(noise)) == pytest.approx(1.5, abs=0.15)


def test_long_window_hrv_is_dormant_on_10_s_and_fires_on_device_length():
    short = extract(synth(), strict=True)
    for name in LONG:
        assert math.isnan(short[name]) and "requires >= 300 s" in short[f"{name}_reason"], name
    device = extract(synth(leads=["II"], duration_s=320), strict=True)
    assert device["sdnn_ms"] == pytest.approx(0.0, abs=2.0)   # regular synthetic rhythm
    assert device["sdnn_ms_reason"] == ""
    # a constant RR series has no spectrum and no entropy: a stated reason, not a number
    assert math.isnan(device["lf_hf_ratio"]) and "variability" in device["lf_hf_ratio_reason"]
