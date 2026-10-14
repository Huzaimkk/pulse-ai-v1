"""
Atrial family tests — written before pulse_ecg/features/atrial.py.

Synthetic truth comes from pulse_ecg.synth: half-sine P waves with exact support,
limb leads Einthoven/Goldberger-consistent (P axis = atan2(2*aVF/sqrt(3), I)),
and an optional biphasic V1 P whose PTFV1 is known in closed form.
"""

import math

import numpy as np
import pytest

from pulse_ecg import config
from pulse_ecg.registry import REGISTRY, NotApplicable
from pulse_ecg.synth import synth

A = config.load()["atrial"]
ATRIAL = ["p_duration_ms", "ptfv1_uv_ms", "p_axis_deg", "iab_class", "mvp_score"]


def extract(rec):
    from pulse_ecg.extract import extract as _extract
    return _extract(rec, strict=True)


def atrial():
    from pulse_ecg.features import atrial as mod
    return mod


# ---------------------------------------------------------- declarations ---
def test_every_atrial_feature_is_registered_and_fully_declared():
    atrial()
    for name in ATRIAL:
        s = REGISTRY[name]
        assert s.family == "atrial"
        assert s.units is not None and s.leads and s.cite and s.desc
        assert s.min_duration_s == A["min_duration_s"]


def test_thresholds_come_from_config():
    atrial()
    assert REGISTRY["p_duration_ms"].normal == tuple(A["p_duration_ms"]["normal"])
    assert str(REGISTRY["p_duration_ms"].abnormal) == ">= 120"
    assert str(REGISTRY["ptfv1_uv_ms"].abnormal) == "<= -4000"
    assert REGISTRY["p_axis_deg"].normal == (0.0, 75.0)
    assert REGISTRY["mvp_score"].normal == tuple(A["mvp"]["normal"])
    assert str(REGISTRY["iab_class"].abnormal) == ">= 1"


def test_primary_and_secondary_status_is_stated():
    mod = atrial()
    assert "PRIMARY ENDPOINT" in mod.p_duration_ms.__doc__
    assert "SECONDARY" in mod.iab_class.__doc__ and "DESCRIPTIVE" in mod.iab_class.__doc__
    assert REGISTRY["ptfv1_uv_ms"].leads == ("V1",)


# ------------------------------------------------------------ P duration ---
@pytest.mark.parametrize("p_ms, pr_ms", [(90, 160), (110, 160), (130, 180), (150, 200)])
def test_p_duration_recovers_known_width(p_ms, pr_ms):
    out = extract(synth(p_duration_ms=p_ms, pr_ms=pr_ms))
    assert out["p_duration_ms_reason"] == ""
    assert abs(out["p_duration_ms"] - p_ms) < 8


def test_p_duration_is_available_on_device_lead():
    out = extract(synth(leads=["II"]))
    assert out["p_duration_ms_reason"] == ""
    assert abs(out["p_duration_ms"] - 100) < 8


# ----------------------------------------------------------------- PTFV1 ---
def test_ptfv1_matches_closed_form_for_biphasic_v1():
    rec, tr = synth(biphasic_v1=True, return_truth=True)
    out = extract(rec)
    assert out["ptfv1_uv_ms"] < 0
    assert abs(out["ptfv1_uv_ms"] - tr["ptfv1_uv_ms"]) <= 0.25 * abs(tr["ptfv1_uv_ms"])


def test_ptfv1_deeper_terminal_part_is_more_negative():
    shallow = extract(synth(biphasic_v1=True, v1_terminal_amp_mv=0.06))["ptfv1_uv_ms"]
    deep = extract(synth(biphasic_v1=True, v1_terminal_amp_mv=0.15))["ptfv1_uv_ms"]
    assert deep < shallow < 0


def test_ptfv1_is_zero_for_upright_v1_p():
    assert extract(synth())["ptfv1_uv_ms"] == 0


def test_terminal_negative_on_array():
    fs = 1000
    x = np.r_[np.sin(np.pi * np.arange(50) / 50) * 0.1,
              -np.sin(np.pi * np.arange(40) / 40) * 0.08]
    dur_ms, depth_uv = atrial().terminal_negative(x, fs, A["min_deflection_mv"])
    assert dur_ms == pytest.approx(40, abs=2)
    assert depth_uv == pytest.approx(80, rel=0.02)
    assert atrial().terminal_negative(np.abs(x), fs, A["min_deflection_mv"]) == (0.0, 0.0)


def test_ptfv1_raises_not_applicable_without_v1():
    from pulse_ecg.extract import context
    ctx = context(synth(leads=["II"]))
    with pytest.raises(NotApplicable, match="V1"):
        atrial().ptfv1_uv_ms(ctx)


# ---------------------------------------------------------------- P axis ---
@pytest.mark.parametrize("axis", [-20, 10, 50, 70, 100])
def test_p_axis_recovers_synthetic_axis(axis):
    out = extract(synth(p_axis_deg=axis))
    assert abs(out["p_axis_deg"] - axis) < 5


def test_p_axis_from_areas_uses_hexaxial_geometry():
    th = math.radians(60)
    i, avf = math.cos(th), math.sqrt(3) / 2 * math.sin(th)
    assert atrial().p_axis_from_areas(i, avf) == pytest.approx(60)


# ------------------------------------------------------------------- IAB ---
def test_biphasic_detection():
    t = np.arange(100)
    pos_neg = np.r_[np.sin(np.pi * t[:50] / 50) * 0.1, -np.sin(np.pi * t[:50] / 50) * 0.1]
    m = A["min_deflection_mv"]
    assert atrial().biphasic_pos_neg(pos_neg, m)
    assert not atrial().biphasic_pos_neg(-pos_neg, m)          # -/+ is not +/-
    assert not atrial().biphasic_pos_neg(np.abs(pos_neg), m)    # monophasic
    assert not atrial().biphasic_pos_neg(np.r_[pos_neg[:50], -0.01 * pos_neg[:50]], m)


def test_iab_classification_rule():
    f = atrial().classify_iab
    assert f(110, 3) == 0          # short P, whatever the morphology
    assert f(125, 0) == 1          # partial
    assert f(125, A["iab"]["advanced_min_leads"]) == 2


def test_iab_class_on_synthetic_p():
    assert extract(synth(p_duration_ms=100))["iab_class"] == 0
    assert extract(synth(p_duration_ms=135, pr_ms=180))["iab_class"] == 1


def test_iab_needs_inferior_leads():
    out = extract(synth(leads=["II"]))
    assert math.isnan(out["iab_class"]) and "requires lead III" in out["iab_class_reason"]


# ------------------------------------------------------------------- MVP ---
def test_mvp_points():
    f = atrial().mvp_points
    assert f(0, 0.25, 100) == 0
    assert f(0, 0.15, 100) == 1
    assert f(0, 0.05, 100) == 2
    assert f(0, 0.25, 130) == 1
    assert f(0, 0.25, 145) == 2
    assert f(2, 0.05, 150) == 6
    assert f(0, 0.2, 120) == 2     # boundaries fall in the middle bands
    assert f(0, 0.1, 140) == 2


P155_XFAIL = ("clean half-sine P >= 155 ms is not delineated at p_span_ms 100; measured "
              "0/981 on real Design A records, see paper/notes/p_missingness.md")


@pytest.mark.parametrize("kw, expect", [
    (dict(p_amp_mv=0.4), 0),                                   # I ~0.26 mV, P 100
    (dict(p_amp_mv=0.23), 1),                                  # I ~0.15 mV
    (dict(p_amp_mv=0.1), 2),                                   # I ~0.06 mV
    pytest.param(dict(p_amp_mv=0.1, p_duration_ms=155, pr_ms=200), 5,  # V 2 + P 2 + partial IAB 1
                 marks=pytest.mark.xfail(strict=True, reason=P155_XFAIL)),
    pytest.param(dict(p_axis_deg=70, p_duration_ms=155, pr_ms=200), 5,  # I ~0.05 mV at default amplitude
                 marks=pytest.mark.xfail(strict=True, reason=P155_XFAIL)),
])
def test_mvp_score_on_synthetic_ecg(kw, expect):
    assert extract(synth(**kw))["mvp_score"] == expect


# ---------------------------------------------------------- device path ---
def test_device_recording_never_crashes_and_states_reasons():
    out = extract(synth(leads=["II"]))
    for name in ["ptfv1_uv_ms", "p_axis_deg", "iab_class", "mvp_score"]:
        assert math.isnan(out[name]), name
        assert out[f"{name}_reason"].startswith("requires lead"), name
    assert "V1" in out["ptfv1_uv_ms_reason"]
