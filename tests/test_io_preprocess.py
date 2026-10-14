"""Readers and filtering. Dataset tests skip when the data is not downloaded."""

import numpy as np
import pytest
from types import SimpleNamespace

from pulse_ecg import config, io
from pulse_ecg.preprocess import filter_signal, preprocess
from pulse_ecg.synth import synth


def _ann(symbols, samples):
    return SimpleNamespace(symbol=symbols, sample=np.asarray(samples))


def test_parse_waves_triplets():
    w = io.parse_waves(_ann(list("(p)(N)(t)"), range(9)), "II")
    assert list(w.wave) == ["p", "qrs", "t"]
    assert w.onset.tolist() == [0, 3, 6] and w.offset.tolist() == [2, 5, 8]


def test_parse_waves_missing_onset_is_nan_not_borrowed():
    # QTDB style: T wave marked 't' ')' with no '('
    w = io.parse_waves(_ann(["(", "N", ")", "t", ")"], [0, 1, 2, 3, 4]), "all")
    t = w[w.wave == "t"].iloc[0]
    assert np.isnan(t.onset) and t.offset == 4


def test_filter_removes_dc_and_50hz_keeps_shape():
    rec = synth()
    t = np.arange(rec.n_samples) / rec.fs
    dirty = rec.signal + 0.5 + 0.2 * np.sin(2 * np.pi * 50 * t)[:, None]
    clean = filter_signal(dirty, rec.fs)
    ref = filter_signal(rec.signal, rec.fs)
    assert np.abs(clean - ref)[500:-500].max() < 0.01


def test_preprocess_keeps_recording_identity():
    rec = synth(leads=["II"])
    out = preprocess(rec)
    assert out.leads == rec.leads and out.fs == rec.fs and out.meta["filtered"]


def test_ptbxl_path_refuses_records100():
    with pytest.raises(ValueError, match="Rule 1"):
        io.ptbxl_path({"filename_hr": "records100/00000/00001_lr"})


needs = lambda key, f: pytest.mark.skipif(  # noqa: E731
    not (config.path(key) / f).exists(), reason=f"{key} not downloaded")


@needs("ludb", "RECORDS")
def test_read_ludb_is_12_lead_500hz_with_waves():
    rec, waves = io.read_ludb("1")
    assert rec.fs == 500 and rec.n_leads == 12 and rec.duration_s == 10
    assert set(waves.lead) == set(rec.leads)
    assert {"p", "qrs", "t"} <= set(waves.wave)


@needs("qtdb", "RECORDS")
def test_read_qtdb_is_250hz():
    rec, waves = io.read_qtdb("sel100")
    assert rec.fs == 250 and rec.n_leads == 2 and len(waves) > 0


@needs("ptbxl", "scp_statements.csv")
def test_scp_statements_have_the_cohort_codes():
    scp = io.scp_statements()
    assert {"LAO/LAE", "SR", "AFIB"} <= set(scp.index)


@needs("ptbxl", "records500/00000/00001_hr.dat")
def test_read_ptbxl_first_record():
    db = io.ptbxl_database()
    rec = io.read_ptbxl(db.loc[1])
    assert rec.fs == 500 and rec.n_samples == 5000 and rec.n_leads == 12
    assert rec.meta["patient_id"] == db.loc[1, "patient_id"]
