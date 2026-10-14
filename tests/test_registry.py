import math

import numpy as np
import pandas as pd
import pytest

from pulse_ecg.recording import STANDARD_LEADS, Recording
from pulse_ecg.registry import Context, Registry, Threshold


def rec(leads=STANDARD_LEADS, seconds=10.0, fs=500):
    n = int(seconds * fs)
    return Recording.from_array(np.zeros((n, len(leads))), fs, leads, "synthetic", "r")


@pytest.fixture
def reg():
    r = Registry()

    @r.feature(name="p_duration_ms", family="atrial", units="ms", leads=["vm"],
               requires=["p_onset", "p_offset"], normal=(80, 120), abnormal=(">=", 120),
               cite="martinezselles2020bayes", desc="P duration")
    def _p(ctx):
        return (ctx.p_offset - ctx.p_onset) / ctx.fs * 1000

    @r.feature(name="ptfv1_uv_ms", family="atrial", units="µV·ms", leads=["V1"],
               abnormal=("<=", -4000), cite="x", desc="PTFV1")
    def _ptf(ctx):
        return -5000.0

    @r.feature(name="sdnn_ms", family="rate", units="ms", leads=["any"],
               min_duration_s=300, desc="SDNN over 5 min")
    def _sdnn(ctx):
        return 40.0

    return r


def test_computes_value_and_empty_reason(reg):
    out = reg.compute(rec(), Context(rec(), {"p_onset": 100, "p_offset": 155}), strict=True)
    assert out["p_duration_ms"] == pytest.approx(110.0)
    assert out["p_duration_ms_reason"] == ""


def test_missing_lead_gives_nan_with_reason(reg):
    out = reg.compute(rec(leads=["II"]))
    assert math.isnan(out["ptfv1_uv_ms"])
    assert out["ptfv1_uv_ms_reason"] == "requires lead V1, not present"


def test_short_recording_gates_long_window_features(reg):
    out = reg.compute(rec())                       # 10 s, like PTB-XL
    assert math.isnan(out["sdnn_ms"])
    assert "requires >= 300 s" in out["sdnn_ms_reason"]
    long_out = reg.compute(rec(leads=["II"], seconds=600, fs=100))   # device-length
    assert long_out["sdnn_ms"] == 40.0


def test_failed_delineation_gives_nan_with_reason(reg):
    out = reg.compute(rec(), Context(rec(), {"p_onset": 100, "p_offset": float("nan")}))
    assert math.isnan(out["p_duration_ms"])
    assert "p_offset" in out["p_duration_ms_reason"]


def test_exception_is_recorded_not_raised_unless_strict():
    r = Registry()

    @r.feature(name="boom", family="quality", units="", leads=["any"], desc="raises")
    def _boom(ctx):
        raise ZeroDivisionError("no beats")

    out = r.compute(rec())
    assert math.isnan(out["boom"]) and out["boom_reason"] == "error: ZeroDivisionError: no beats"
    with pytest.raises(ZeroDivisionError):
        r.compute(rec(), strict=True)


def test_declaration_is_validated():
    r = Registry()
    ok = dict(family="atrial", units="ms", leads=["II"], desc="d")
    r.feature(name="a", **ok)(lambda c: 1)
    with pytest.raises(ValueError, match="already registered"):
        r.feature(name="a", **ok)(lambda c: 1)
    with pytest.raises(ValueError, match="unknown lead"):
        r.feature(name="b", family="atrial", units="ms", leads=["V7"], desc="d")
    with pytest.raises(ValueError, match="family"):
        r.feature(name="c", family="cardiac", units="ms", leads=["II"], desc="d")
    with pytest.raises(ValueError, match="citation"):
        r.feature(name="d", normal=(1, 2), **ok)
    with pytest.raises(ValueError, match="lo < hi"):
        r.feature(name="e", normal=(2, 1), cite="k", **ok)
    with pytest.raises(ValueError, match="snake_case"):
        r.feature(name="PDur", **ok)
    with pytest.raises(ValueError, match="operator"):
        Threshold("=>", 1)


def test_abnormal_threshold_evaluates_and_prints(reg):
    assert reg.is_abnormal("p_duration_ms", 125.0) is True
    assert reg.is_abnormal("p_duration_ms", 100.0) is False
    assert reg.is_abnormal("p_duration_ms", float("nan")) is None
    assert str(reg["ptfv1_uv_ms"].abnormal) == "<= -4000"


def test_capability_map_for_device(reg):
    cap = reg.capability(["II"], duration_s=3600).set_index("feature")
    assert cap.loc["p_duration_ms", "available"]          # VM from one lead: partial
    assert not cap.loc["ptfv1_uv_ms", "available"]
    assert cap.loc["sdnn_ms", "available"]                # long device recording
    ptb = reg.capability(STANDARD_LEADS, duration_s=10).set_index("feature")
    assert not ptb.loc["sdnn_ms", "available"]


def test_table_is_derived_from_declarations(reg):
    t = reg.table().set_index("feature")
    assert t.loc["p_duration_ms", "normal"] == "80–120"
    assert t.loc["p_duration_ms", "abnormal"] == ">= 120"
    assert t.loc["p_duration_ms", "cite"] == "martinezselles2020bayes"
    assert t.loc["sdnn_ms", "abnormal"] == ""


def test_sanity_reports_fraction_in_normal_range(reg):
    df = pd.DataFrame({"p_duration_ms": [90, 100, 130, np.nan]})
    s = reg.sanity(df).set_index("feature")
    assert s.loc["p_duration_ms", "n"] == 3
    assert s.loc["p_duration_ms", "frac_in_normal"] == pytest.approx(2 / 3)


def test_not_applicable_is_a_reason_not_an_error():
    from pulse_ecg.registry import NotApplicable
    r = Registry()

    @r.feature(name="fqrs", family="ventricular", units="", leads=["any"], desc="narrow QRS only")
    def _f(ctx):
        raise NotApplicable("QRS >= 120 ms")

    for strict in (False, True):                    # declining is not a bug: never raised
        out = r.compute(rec(), strict=strict)
        assert math.isnan(out["fqrs"]) and out["fqrs_reason"] == "QRS >= 120 ms"
