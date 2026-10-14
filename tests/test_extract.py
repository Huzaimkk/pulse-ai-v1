"""The one entry point: Recording -> one flat row (build order §7, before any feature)."""

import math

import numpy as np
import pytest

import pulse_ecg.features  # noqa: F401  (registers features)
from pulse_ecg.delineate import delineate
from pulse_ecg.extract import extract
from pulse_ecg.preprocess import preprocess
from pulse_ecg.registry import REGISTRY
from pulse_ecg.report import prf
from pulse_ecg.synth import synth


def test_row_has_identity_and_a_reason_for_every_feature():
    rec = synth()
    out = extract(rec, strict=True)
    assert out["record_id"] == rec.record_id and out["source"] == "synthetic"
    assert out["n_leads"] == 12 and out["duration_s"] == pytest.approx(10.0)
    assert out["delineation_reason"] == ""
    for spec in REGISTRY:
        assert spec.name in out and spec.reason_key in out
        nan = isinstance(out[spec.name], float) and math.isnan(out[spec.name])
        assert nan == bool(out[spec.reason_key]), spec.name     # NaN <=> stated reason


def test_failed_delineation_gives_every_beat_feature_a_reason():
    out = extract(synth(duration_s=1.2))
    assert "complete beats" in out["delineation_reason"]
    assert math.isnan(out["qrs_duration_ms"]) and out["qrs_duration_ms_reason"]


def test_single_lead_marks_are_kept_per_lead():
    rec = synth(biphasic_v1=True)
    d = delineate(preprocess(rec))
    assert set(d.per_lead) == set(rec.leads)
    v1 = d.per_lead["V1"]
    assert v1["p_offset"] > v1["p_onset"]             # median-beat samples, V1 alone


def test_prf_matches_definition_and_handles_empty_classes():
    y = np.array([1, 1, 1, 0, 0, 0, 1, 0])
    yhat = np.array([1, 1, 0, 1, 0, 0, 1, 0])
    p, r, f1 = prf(y, yhat)                           # tp 3, fp 1, fn 1
    assert (p, r) == (0.75, 0.75) and f1 == pytest.approx(0.75)
    p, r, f1 = prf(np.array([1, 0]), np.array([0, 0]))  # nothing predicted positive
    assert math.isnan(p) and r == 0.0 and math.isnan(f1)
