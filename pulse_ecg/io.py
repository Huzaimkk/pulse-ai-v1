"""
Readers. Every source becomes a Recording; reference annotations (LUDB, QTDB)
become a tidy wave table. Nothing downstream imports wfdb.

    read_ptbxl(row)          -> Recording                 (records500 only, Rule 1)
    read_ludb(n)             -> Recording, waves          (Gate 2, primary)
    read_qtdb(name)          -> Recording, waves          (Gate 2, secondary)
    ptbxl_database()         -> DataFrame, scp_codes parsed to dicts
    scp_statements()         -> DataFrame indexed by SCP code

A wave table has one row per annotated wave per lead:
    lead, wave ("p" | "qrs" | "t"), onset, peak, offset     (samples; NaN if absent)
"""

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import wfdb

from . import config
from .recording import Recording, canonical_lead

WAVES = {"p": "p", "N": "qrs", "t": "t"}
WAVE_COLUMNS = ["lead", "wave", "onset", "peak", "offset"]


# ------------------------------------------------------------------ PTB-XL ---
def ptbxl_database() -> pd.DataFrame:
    root = config.path("ptbxl")
    db = pd.read_csv(root / "ptbxl_database.csv", index_col="ecg_id")
    db["scp_codes"] = db["scp_codes"].apply(ast.literal_eval)
    db["patient_id"] = db["patient_id"].astype(int)
    return db


def scp_statements() -> pd.DataFrame:
    return pd.read_csv(config.path("ptbxl") / "scp_statements.csv", index_col=0)


def ptbxl_path(row) -> Path:
    """Path stem of a record's 500 Hz file. Rule 1: never records100."""
    rel = row["filename_hr"]
    if not str(rel).startswith(config.load()["ptbxl"]["records"] + "/"):
        raise ValueError(f"{rel}: not a records500 file (Rule 1)")
    return config.path("ptbxl") / rel


def read_ptbxl(row, ecg_id=None) -> Recording:
    """`row` is a ptbxl_database() row (its .name is the ecg_id)."""
    ecg_id = row.name if ecg_id is None else ecg_id
    sig, fields = wfdb.rdsamp(str(ptbxl_path(row)))
    _check_units(fields, f"ptbxl {ecg_id}")
    meta = {k: row[k] for k in ("patient_id", "age", "sex") if k in row}
    return Recording.from_array(sig, fields["fs"], fields["sig_name"], "ptbxl",
                                str(ecg_id), meta)


# -------------------------------------------------------------------- LUDB ---
LUDB_LEAD_EXT = {"I": "i", "II": "ii", "III": "iii", "aVR": "avr", "aVL": "avl",
                 "aVF": "avf", "V1": "v1", "V2": "v2", "V3": "v3", "V4": "v4",
                 "V5": "v5", "V6": "v6"}


def ludb_records() -> list[str]:
    return (config.path("ludb") / "RECORDS").read_text().split()


def read_ludb(record) -> tuple[Recording, pd.DataFrame]:
    """record: '1'..'200' or 'data/1'. Waves are annotated per lead."""
    stem = str(record).split("/")[-1]
    base = config.path("ludb") / "data" / stem
    sig, fields = wfdb.rdsamp(str(base))
    _check_units(fields, f"ludb {stem}")
    rec = Recording.from_array(sig, fields["fs"], fields["sig_name"], "ludb", stem,
                               {"comments": fields.get("comments", [])})
    tables = [parse_waves(wfdb.rdann(str(base), LUDB_LEAD_EXT[lead]), lead)
              for lead in rec.leads]
    return rec, pd.concat(tables, ignore_index=True)


# -------------------------------------------------------------------- QTDB ---
def qtdb_records() -> list[str]:
    return (config.path("qtdb") / "RECORDS").read_text().split()


def read_qtdb(record, annotator="q1c") -> tuple[Recording, pd.DataFrame]:
    """Two-lead Holter excerpt; `q1c` is cardiologist 1's manual annotation.
    Annotations are not per lead, so the wave table's lead is 'all'."""
    base = config.path("qtdb") / str(record)
    sig, fields = wfdb.rdsamp(str(base))
    _check_units(fields, f"qtdb {record}")
    leads = _unique_names(fields["sig_name"])
    rec = Recording.from_array(sig, fields["fs"], leads, "qtdb", str(record))
    return rec, parse_waves(wfdb.rdann(str(base), annotator), "all")


# ----------------------------------------------------------------- helpers ---
def parse_waves(ann, lead: str) -> pd.DataFrame:
    """WFDB '(' peak ')' triplets -> one row per wave.

    A boundary belongs to a peak only if no other peak lies between them, so a
    missing '(' (QTDB often marks T waves as 't' ')') gives NaN, not the onset
    of a neighbouring wave."""
    sym, smp = list(ann.symbol), np.asarray(ann.sample)
    rows = []
    for i, s in enumerate(sym):
        if s not in WAVES:
            continue
        onset = offset = np.nan
        if i > 0 and sym[i - 1] == "(":
            onset = smp[i - 1]
        if i + 1 < len(sym) and sym[i + 1] == ")":
            offset = smp[i + 1]
        rows.append((lead, WAVES[s], onset, smp[i], offset))
    return pd.DataFrame(rows, columns=WAVE_COLUMNS).astype(
        {"onset": float, "peak": float, "offset": float})


def _check_units(fields, name):
    units = {str(u).lower() for u in fields.get("units", [])}
    if units and units != {"mv"}:
        raise ValueError(f"{name}: expected mV, got units {sorted(units)}")


def _unique_names(names):
    out, seen = [], {}
    for n in names:
        n = canonical_lead(n)
        seen[n] = seen.get(n, 0) + 1
        out.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    return out
