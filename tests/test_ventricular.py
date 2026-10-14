"""Ventricular family: known answers on synthetic ECGs.

Interval tests compare two synthetic widths, so the fixture's constant boundary
bias (tests/test_delineate.py) cancels and what is tested is the feature, not the
delineator a second time."""

import math

import numpy as np
import pytest

from pulse_ecg.extract import context, extract
from pulse_ecg.features import ventricular
from pulse_ecg.recording import Recording
from pulse_ecg.synth import synth

FS = 500


def run(**kw):
    return extract(synth(**kw), strict=True)


def test_qrs_duration_tracks_the_synthetic_width():
    a, b = run(qrs_duration_ms=90), run(qrs_duration_ms=115)
    assert b["qrs_duration_ms"] - a["qrs_duration_ms"] == pytest.approx(25, abs=8)


def test_qt_tracks_the_synthetic_interval():
    a, b = run(qt_ms=380), run(qt_ms=440)
    assert b["qt_ms"] - a["qt_ms"] == pytest.approx(60, abs=8)


@pytest.mark.parametrize("hr", [60, 90])
def test_qtc_formulas(hr):
    out = run(hr=hr, qt_ms=380 if hr == 90 else 400)
    rr_s = np.mean(out["rr_intervals_ms"]) / 1000
    assert out["qtc_bazett_ms"] == pytest.approx(out["qt_ms"] / math.sqrt(rr_s))
    assert out["qtc_fridericia_ms"] == pytest.approx(out["qt_ms"] / rr_s ** (1 / 3))
    if hr == 60:                                      # RR = 1 s: every correction is QT
        assert out["qtc_bazett_ms"] == pytest.approx(out["qt_ms"], abs=1)


def test_tpeak_tend_matches_the_synthetic_t_wave():
    # T is a half-sine ending at QRS onset + QT, width (QT - QRS) / 2, so its peak
    # is (QT - QRS) / 4 before its end: 77.5 ms for QRS 90, QT 400.
    out = run()
    assert out["tpe_ms"] == pytest.approx(77.5, abs=12)
    assert 0.15 <= out["tpe_qt_ratio"] <= 0.25
    assert out["tpe_qt_ratio"] == pytest.approx(out["tpe_ms"] / out["qt_ms"])


# ------------------------------------------------------------------ fQRS ---
def _complex(*waves, n=60):
    """Chain of half-sine deflections (amplitude, width in samples) -> 1-D QRS."""
    x = [np.zeros(5)]
    for amp, w in waves:
        x.append(amp * np.sin(np.pi * np.arange(w) / w))
    x.append(np.zeros(5))
    return np.concatenate(x)


@pytest.mark.parametrize("waves, frag", [
    ([(-0.1, 6), (1.0, 20), (-0.3, 12)], False),              # qRs
    ([(1.0, 20), (-0.5, 14)], False),                         # RS
    ([(-0.8, 30)], False),                                    # QS
    ([(0.3, 10), (-0.8, 16), (0.5, 12)], True),               # rSR'
    ([(1.0, 20), (-0.5, 8), (0.15, 5), (-0.4, 8)], True),     # notched S nadir
    ([(-0.6, 12), (0.1, 5), (-0.6, 12)], True),               # notched QS (W)
])
def test_fragmentation_of_a_single_complex(waves, frag):
    assert ventricular.fragmented(_complex(*waves), min_prominence_mv=0.05) is frag


def _notched(leads, **kw):
    """Synthetic ECG with an extra deflection (R') at the end of the QRS in `leads`."""
    rec, truth = synth(return_truth=True, **kw)
    sig = rec.signal.copy()
    w = int(0.012 * FS)
    bump = 0.25 * np.sin(np.pi * np.arange(w) / w)
    for q_off in truth["qrs_offset"]:
        for lead in leads:
            j = rec.leads.index(lead)
            sig[q_off - 2 * w:q_off - w, j] += bump
    return Recording.from_array(sig, FS, rec.leads, "synthetic", "notched")


def _fragmented_leads(rec):
    ctx = context(rec)
    beat = ventricular._wideband_beat(ctx)
    lo, hi = int(round(ctx.points["qrs_onset"])), int(round(ctx.points["qrs_offset"]))
    return {l for l in rec.leads if ventricular.fragmented(beat[lo:hi + 1, rec.leads.index(l)])}


def test_v1_rsr_prime_counts_as_fragmented():
    # The fixture's V1 is r (+) S (-) r' (+): an rSr' pattern, which the Das notch
    # criterion counts as fragmented. Expected behaviour, not a bug: rSr' in V1
    # occurs in healthy subjects, so anterior-territory fQRS on real data inherits
    # this false positive (paper/notes/fqrs_v1_rsr_prime.md, Limitations).
    assert "V1" in _fragmented_leads(synth())


def test_fqrs_needs_two_leads_of_one_territory():
    assert run()["fqrs_present"] == 0
    assert extract(_notched(["II", "III", "aVF"]), strict=True)["fqrs_present"] == 1


def test_one_fragmented_lead_alone_in_its_territory_is_not_fqrs():
    # True negative control: one notch, in II; no other inferior lead is fragmented,
    # so no territory reaches min_leads. II itself must be detected, so the 0 comes
    # from the territory rule and not from a missed notch.
    rec = _notched(["II"])
    frag = _fragmented_leads(rec)
    assert "II" in frag and not frag & {"III", "aVF"}
    assert extract(rec, strict=True)["fqrs_present"] == 0


def test_fqrs_does_not_apply_to_a_wide_qrs():
    out = run(qrs_duration_ms=150, qt_ms=460)
    assert math.isnan(out["fqrs_present"]) and "120" in out["fqrs_present_reason"]


# ------------------------------------------------------------ ST at J+60 ---
def _st_plateau(lead, mv, **kw):
    rec, truth = synth(return_truth=True, **kw)
    sig = rec.signal.copy()
    j = rec.leads.index(lead)
    ramp = int(0.02 * FS)
    for q_off in truth["qrs_offset"]:
        seg = np.full(int(0.135 * FS), mv)                 # J to 20 ms before T onset
        seg[:ramp] = np.linspace(0, mv, ramp)
        seg[-ramp:] = np.linspace(mv, 0, ramp)
        sig[q_off:q_off + len(seg), j] += seg
    return Recording.from_array(sig, FS, rec.leads, "synthetic", "st")


def test_st_level_at_j60():
    flat = run()
    for lead in ("i", "ii", "v2", "v5"):
        assert abs(flat[f"st_j60_{lead}_mv"]) < 0.02, lead
    out = extract(_st_plateau("V2", 0.15), strict=True)
    assert out["st_j60_v2_mv"] == pytest.approx(0.15, abs=0.03)
    assert abs(out["st_j60_v5_mv"]) < 0.02


def test_device_gets_only_the_leads_it_has():
    out = run(leads=["II"])
    assert not math.isnan(out["qrs_duration_ms"]) and not math.isnan(out["st_j60_ii_mv"])
    assert math.isnan(out["st_j60_v2_mv"]) and "requires lead V2" in out["st_j60_v2_mv_reason"]
    assert math.isnan(out["fqrs_present"]) and "requires lead" in out["fqrs_present_reason"]
