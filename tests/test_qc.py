import copy
import inspect
import math

import pandas as pd
import pytest

from pulse_ecg import config, qc
from pulse_ecg.registry import REGISTRY

NAN = math.nan


@pytest.fixture(scope="module")
def cfg():
    return config.load()


@pytest.fixture(scope="module")
def rules(cfg):
    return qc.build_rules(cfg)


def good_row(**kw):
    row = {"delineation_reason": "", "mean_hr_bpm": 70.0, "qrs_duration_ms": 95.0,
           "qt_ms": 400.0, "p_duration_ms": 110.0, "p_offset": 150.0, "qrs_onset": 170.0}
    row.update(kw)
    return row


# ------------------------------------------------------------ derivation ---
def test_ranges_come_from_config_qc_block(cfg, rules):
    got = {r.key: (r.lo, r.hi) for r in rules.ranges}
    assert got == {k: tuple(map(float, v)) for k, v in cfg["qc"].items()}
    assert got == {"hr": (40, 200), "qrs_ms": (50, 200), "qt_ms": (250, 600), "p_ms": (50, 200)}


def test_rules_bound_registered_features_with_registry_units(rules):
    for r in rules.ranges:
        assert r.feature == qc.QC_FEATURES[r.key]
        assert r.units == REGISTRY[r.feature].units


def test_ranges_contain_registry_normal_range_and_abnormal_threshold(rules):
    for r in rules.ranges:
        spec = REGISTRY[r.feature]
        if spec.normal:
            assert r.lo <= spec.normal[0] and spec.normal[1] <= r.hi
        if spec.abnormal:
            assert r.lo <= spec.abnormal.value <= r.hi


def test_range_that_cuts_into_normal_range_is_rejected(cfg):
    bad = copy.deepcopy(cfg)
    bad["qc"]["p_ms"] = [90, 200]                  # p_duration normal is 80-120
    with pytest.raises(ValueError, match="normal range"):
        qc.build_rules(bad)


def test_range_that_excludes_abnormal_threshold_is_rejected(cfg):
    bad = copy.deepcopy(cfg)
    bad["qc"]["qrs_ms"] = [50, 115]                # qrs abnormal is >= 120
    with pytest.raises(ValueError):
        qc.build_rules(bad)


def test_unmapped_qc_key_is_rejected(cfg):
    bad = copy.deepcopy(cfg)
    bad["qc"]["foo_ms"] = [0, 1]
    with pytest.raises(ValueError, match="foo_ms"):
        qc.build_rules(bad)


def test_config_change_moves_the_range(cfg):
    alt = copy.deepcopy(cfg)
    alt["qc"]["hr"] = [30, 220]
    r = {x.key: x for x in qc.build_rules(alt).ranges}["hr"]
    assert (r.lo, r.hi) == (30, 220)


def test_required_features_are_range_features_plus_primary_endpoint(cfg, rules):
    assert set(rules.required) == set(qc.QC_FEATURES.values()) | {cfg["analysis"]["primary"]["feature"]}


# ------------------------------------------------------------ evaluation ---
def test_good_record_passes(rules):
    res = qc.evaluate(good_row(), rules)
    assert res.passed and res.failed == []


@pytest.mark.parametrize("feature,value,rule", [
    ("mean_hr_bpm", 39.9, "range:hr"), ("mean_hr_bpm", 200.1, "range:hr"),
    ("qrs_duration_ms", 49.0, "range:qrs_ms"), ("qrs_duration_ms", 201.0, "range:qrs_ms"),
    ("qt_ms", 249.0, "range:qt_ms"), ("qt_ms", 601.0, "range:qt_ms"),
    ("p_duration_ms", 49.0, "range:p_ms"), ("p_duration_ms", 201.0, "range:p_ms"),
])
def test_out_of_range_fails_named_rule(rules, feature, value, rule):
    res = qc.evaluate(good_row(**{feature: value}), rules)
    assert not res.passed and res.rules == [rule]
    assert feature in res.failed[0][1]


@pytest.mark.parametrize("feature,value", [
    ("mean_hr_bpm", 40.0), ("mean_hr_bpm", 200.0), ("qrs_duration_ms", 50.0),
    ("qrs_duration_ms", 200.0), ("qt_ms", 250.0), ("qt_ms", 600.0),
    ("p_duration_ms", 50.0), ("p_duration_ms", 200.0)])
def test_range_bounds_are_inclusive(rules, feature, value):
    assert qc.evaluate(good_row(**{feature: value}), rules).passed


def test_qrs_onset_before_p_offset_fails(rules):
    res = qc.evaluate(good_row(p_offset=160.0, qrs_onset=150.0), rules)
    assert res.rules == [qc.OVERLAP_RULE]
    assert "150" in res.failed[0][1] and "160" in res.failed[0][1]


def test_qrs_onset_equal_to_p_offset_passes(rules):
    assert qc.evaluate(good_row(p_offset=160.0, qrs_onset=160.0), rules).passed


def test_overlap_fails_even_when_every_feature_is_in_range(rules):
    # the rule exists because overlap corrupts features without pushing them out of range
    res = qc.evaluate(good_row(p_offset=160.0, qrs_onset=159.5), rules)
    assert not res.passed


def test_nan_required_feature_fails_with_reason_carried(rules):
    row = good_row(p_duration_ms=NAN, p_duration_ms_reason="requires p_onset, not delineated")
    res = qc.evaluate(row, rules)
    assert res.failed == [("missing:p_duration_ms", "requires p_onset, not delineated")]


def test_nan_without_reason_still_fails(rules):
    res = qc.evaluate(good_row(qt_ms=NAN), rules)
    assert res.failed == [("missing:qt_ms", "value is NaN")]


@pytest.mark.parametrize("point", qc.FIDUCIALS)
def test_missing_fiducial_fails(rules, point):
    res = qc.evaluate(good_row(**{point: NAN}), rules)
    assert res.rules == [f"missing:{point}"]


def test_failed_delineation_fails_with_reason(rules):
    res = qc.evaluate(good_row(delineation_reason="wavedet found no QRS boundaries"), rules)
    assert ("delineation", "wavedet found no QRS boundaries") in res.failed


def test_every_failed_rule_is_listed(rules):
    row = good_row(mean_hr_bpm=250.0, qt_ms=NAN, p_offset=180.0, qrs_onset=170.0)
    assert set(qc.evaluate(row, rules).rules) == {"range:hr", "missing:qt_ms", qc.OVERLAP_RULE}


def test_evaluate_is_blind_to_group(rules):
    assert list(inspect.signature(qc.evaluate).parameters) == ["row", "rules"]
    a = qc.evaluate(good_row(group="case", p_offset=160.0, qrs_onset=150.0), rules)
    b = qc.evaluate(good_row(group="control", p_offset=160.0, qrs_onset=150.0), rules)
    assert a.failed == b.failed


# ------------------------------------------------------- apply and flow ---
def _frame():
    return pd.DataFrame([
        good_row(ecg_id=1, group="case"),
        good_row(ecg_id=2, group="case", p_offset=160.0, qrs_onset=150.0, qt_ms=700.0),
        good_row(ecg_id=3, group="control", p_duration_ms=NAN),
        good_row(ecg_id=4, group="control"),
        good_row(ecg_id=5, group="control", qrs_duration_ms=220.0),
    ])


def test_apply_keeps_every_record_in_order(rules):
    df = _frame()
    out = qc.apply(df, rules)
    assert len(out) == len(df) and out["ecg_id"].tolist() == df["ecg_id"].tolist()
    assert out["qc_pass"].tolist() == [True, False, False, True, False]
    assert out.loc[1, "qc_failed"] == f"{qc.OVERLAP_RULE};range:qt_ms"
    assert out.loc[1, "qc_first"] == qc.OVERLAP_RULE
    assert (out.loc[~out["qc_pass"], "qc_detail"] != "").all()


def test_flow_accounts_for_every_record(rules):
    f = qc.flow(qc.apply(_frame(), rules), rules, by="group")
    for g, n in [("all", 5), ("case", 2), ("control", 3)]:
        s = f[f["group"] == g].set_index("rule")
        per_rule = s.drop(["records", "excluded", "final_n"])
        assert per_rule["first"].sum() == s.loc["excluded", "first"]
        assert s.loc["excluded", "first"] + s.loc["final_n", "first"] == n
    a = f[f["group"] == "all"].set_index("rule")
    assert a.loc["excluded", "any"] == 3 and a.loc["final_n", "any"] == 2
    assert a.loc[qc.OVERLAP_RULE, "any"] == 1 and a.loc["range:qt_ms", "any"] == 1
    assert a.loc["range:qt_ms", "first"] == 0          # attributed to the overlap rule
    assert a.loc["missing:p_duration_ms", "any"] == 1
    assert a.loc["range:qrs_ms", "any"] == 1


def test_flow_rows_cover_every_rule_name(rules):
    f = qc.flow(qc.apply(_frame(), rules), rules)
    assert set(rules.names) <= set(f["rule"])
