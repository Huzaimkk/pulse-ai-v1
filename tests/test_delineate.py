"""Delineation recovers known boundaries on synthetic ECGs (build order §7, step 3).

Tolerances are tighter than the CSE-style limits Gate 2 is judged against,
because synthetic waves have exact support and no noise unless asked."""

import math

import numpy as np
import pytest

from pulse_ecg.delineate import delineate, detect_r_peaks
from pulse_ecg.preprocess import preprocess
from pulse_ecg.synth import synth

TOL_MS = {"p_onset": 6, "p_offset": 6, "qrs_onset": 10, "qrs_offset": 10, "t_offset": 10}
CASES = [dict(), dict(p_duration_ms=110), dict(hr=90, qt_ms=360),
         dict(hr=45, pr_ms=220, p_duration_ms=130, qt_ms=450),
         dict(leads=["II"]), dict(noise_mv=0.02, seed=3)]


def _run(**kw):
    rec, truth = synth(return_truth=True, **kw)
    return rec, truth, delineate(preprocess(rec))


@pytest.mark.parametrize("kw", CASES, ids=lambda kw: ",".join(f"{k}={v}" for k, v in kw.items()) or "default")
def test_r_peaks_found_once_per_beat(kw):
    rec, truth, d = _run(**kw)
    assert len(d.r_peaks) == truth["n_beats"]
    assert np.all(np.abs(d.r_peaks - truth["qrs_onset"]) < 0.1 * rec.fs)


@pytest.mark.parametrize("kw", CASES, ids=lambda kw: ",".join(f"{k}={v}" for k, v in kw.items()) or "default")
def test_boundaries_match_truth(kw):
    rec, truth, d = _run(**kw)
    assert d.ok, d.reason
    got = d.beat_points()
    for name, tol in TOL_MS.items():
        err_ms = (got[name] - truth[name]) * 1000 / rec.fs
        assert np.all(np.abs(err_ms) < tol), (name, err_ms)


def test_p_duration_on_vector_magnitude():
    rec, truth, d = _run(p_duration_ms=110)
    p = d.points
    assert abs((p["p_offset"] - p["p_onset"]) * 1000 / rec.fs - 110) < 8


def test_no_p_wave_is_reported_absent_not_invented():
    rec, truth, d = _run(p_amp_mv=0.0)
    assert d.ok
    assert math.isnan(d.points["p_onset"]) and math.isnan(d.points["p_offset"])
    assert not math.isnan(d.points["qrs_onset"])


def test_too_short_recording_fails_with_reason():
    rec = synth(duration_s=1.2)
    d = delineate(preprocess(rec))
    assert not d.ok and "complete beats" in d.reason


def test_flat_signal_has_no_r_peaks():
    assert len(detect_r_peaks(np.zeros((5000, 12)), 500)) == 0
