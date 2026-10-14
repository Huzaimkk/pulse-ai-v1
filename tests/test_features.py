"""
Feature tests from Code_Architecture_Plan.md §5 — written BEFORE the features.

Each is xfail(strict=True) until its feature exists: while pending it is
reported as XFAIL; once the feature makes it pass, strict mode turns the
unexpected pass into a FAILURE, forcing the marker to be removed so the test
starts guarding for real.
"""

import math

import numpy as np
import pytest

from pulse_ecg.synth import synth

pending = pytest.mark.xfail(strict=True, reason="feature not implemented yet (build order §7)")


def extract(rec):
    from pulse_ecg.extract import extract as _extract
    return _extract(rec, strict=True)



def test_qrs_duration_in_physiological_range():
    assert 70 <= extract(synth())["qrs_duration_ms"] <= 120



def test_f1_never_exceeds_precision_and_recall():
    # Guards the arithmetic error that reached the previous thesis
    # (precision 93.0 %, recall 92.5 %, F1 93.9 %).
    from pulse_ecg.report import prf
    rng = np.random.default_rng(0)
    for _ in range(200):
        y = rng.integers(0, 2, 300)
        yhat = rng.integers(0, 2, 300)
        p, r, f1 = prf(y, yhat)
        assert f1 <= max(p, r) + 1e-12
        assert min(p, r) - 1e-12 <= f1
