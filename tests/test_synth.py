"""The synthetic generator is the ground truth for every other test — check it first."""

import numpy as np
import pytest

from pulse_ecg.synth import synth


def test_fiducials_match_requested_intervals():
    rec, tr = synth(hr=75, p_duration_ms=110, pr_ms=170, qrs_duration_ms=96,
                    qt_ms=390, return_truth=True)
    ms = lambda a, b: (tr[b] - tr[a]) * 1000 / rec.fs  # noqa: E731
    assert np.allclose(ms("p_onset", "p_offset"), 110, atol=2)
    assert np.allclose(ms("p_onset", "qrs_onset"), 170, atol=2)
    assert np.allclose(ms("qrs_onset", "qrs_offset"), 96, atol=2)
    assert np.allclose(ms("qrs_onset", "t_offset"), 390, atol=2)
    assert 10 <= tr["n_beats"] <= 13          # 75 bpm over ~10 s


def test_signal_is_zero_outside_waves():
    rec, tr = synth(return_truth=True)
    ii = rec.lead("II")
    # isoelectric between T offset of beat k and P onset of beat k+1
    for t_off, p_on in zip(tr["t_offset"][:-1], tr["p_onset"][1:]):
        assert np.all(ii[t_off + 1:p_on] == 0)


def test_limb_leads_obey_einthoven_and_goldberger():
    rec = synth()
    I, II = rec.lead("I"), rec.lead("II")
    assert np.allclose(rec.lead("III"), II - I, atol=1e-12)
    assert np.allclose(rec.lead("aVR"), -(I + II) / 2, atol=1e-12)
    assert np.allclose(rec.lead("aVL"), I - II / 2, atol=1e-12)
    assert np.allclose(rec.lead("aVF"), II - I / 2, atol=1e-12)


@pytest.mark.parametrize("axis", [-20.0, 0.0, 50.0, 75.0, 110.0])
def test_p_axis_is_recoverable_from_i_and_avf(axis):
    rec, tr = synth(p_axis_deg=axis, return_truth=True)
    k = (tr["p_onset"][0] + tr["p_offset"][0]) // 2          # P-wave peak
    i, avf = rec.lead("I")[k], rec.lead("aVF")[k]
    got = np.degrees(np.arctan2(2 * avf / np.sqrt(3), i))
    assert got == pytest.approx(axis, abs=1e-6)


def test_biphasic_v1_has_negative_terminal_p():
    rec, tr = synth(biphasic_v1=True, v1_terminal_amp_mv=0.08, return_truth=True)
    v1 = rec.lead("V1")
    p0, p1 = tr["p_onset"][0], tr["p_offset"][0]
    mid = (p0 + p1) // 2
    assert v1[p0:mid].max() > 0 and v1[mid:p1].min() == pytest.approx(-0.08, abs=1e-3)
    assert tr["ptfv1_uv_ms"] == pytest.approx(-80 * 50)       # 80 µV x 50 ms


def test_single_lead_device_emulation():
    rec = synth(leads=["II"])
    assert rec.leads == ("II",) and rec.n_leads == 1


def test_rejects_impossible_timing():
    with pytest.raises(ValueError):
        synth(p_duration_ms=180, pr_ms=160)
