import numpy as np
import pytest

from pulse_ecg.recording import STANDARD_LEADS, Recording, canonical_lead


def make(n=5000, leads=STANDARD_LEADS, fs=500, **kw):
    return Recording.from_array(np.zeros((n, len(leads))), fs, leads, "synthetic", "r1", **kw)


def test_from_array_derives_duration():
    rec = make()
    assert rec.duration_s == 10.0
    assert (rec.n_samples, rec.n_leads) == (5000, 12)


def test_lead_names_are_canonicalised():
    rec = make(leads=["i", "II", "AVF"], n=10)
    assert rec.leads == ("I", "II", "aVF")
    assert canonical_lead("avr") == "aVR"


def test_nonstandard_lead_names_pass_through():
    # QT Database uses names like MLII; they are allowed but match no feature.
    assert make(leads=["MLII"], n=10).leads == ("MLII",)


def test_rejects_column_lead_mismatch():
    with pytest.raises(ValueError, match="columns"):
        Recording.from_array(np.zeros((10, 3)), 500, ["I", "II"], "synthetic", "r")


def test_rejects_duplicate_leads():
    with pytest.raises(ValueError, match="duplicate"):
        make(leads=["II", "ii"], n=10)


def test_rejects_1d_signal():
    with pytest.raises(ValueError, match="2-D"):
        Recording(np.zeros(10), 500, ("II",), 0.02, "synthetic", "r")


def test_rejects_inconsistent_duration():
    with pytest.raises(ValueError, match="duration_s"):
        Recording(np.zeros((5000, 1)), 500, ("II",), 9.0, "synthetic", "r")


def test_rejects_unknown_source_and_bad_fs():
    with pytest.raises(ValueError, match="source"):
        Recording.from_array(np.zeros((10, 1)), 500, ["II"], "kaggle", "r")
    with pytest.raises(ValueError, match="fs"):
        Recording(np.zeros((10, 1)), 0, ("II",), 0.0, "synthetic", "r")


def test_lead_access_and_select():
    sig = np.arange(24, dtype=float).reshape(2, 12)
    rec = Recording.from_array(sig, 500, STANDARD_LEADS, "synthetic", "r")
    assert rec.lead("V1").tolist() == [6.0, 18.0]
    dev = rec.select(["II"])
    assert dev.leads == ("II",) and dev.lead("II").tolist() == [1.0, 13.0]
    assert dev.duration_s == rec.duration_s
    with pytest.raises(KeyError, match="V1"):
        dev.lead("V1")
    with pytest.raises(KeyError):
        dev.select(["V1"])
