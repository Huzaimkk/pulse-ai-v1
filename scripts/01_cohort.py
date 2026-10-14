"""
01_cohort.py — Gate 1: build both cohort designs from PTB-XL and count them.
Paper Table I.

    python scripts/01_cohort.py
    python scripts/01_cohort.py --min-likelihood 100     # sensitivity check

Every analysed record is in sinus rhythm (P waves are needed). Two designs:
  A (primary)    cases: sinus records carrying LAO/LAE
                 controls: sinus records of patients who never carry LAO/LAE
  B (secondary)  cases: sinus records of patients with AFIB in another record
                 controls: sinus records of patients who never carry AFIB
One record per patient (the earliest eligible), so no heart appears twice
(Rule 2). Controls are matched exactly on cohort.match_exact (sex) and on age within cohort.age_tol years.

SCP likelihood 0 means "stated, likelihood not given", not "absent", so the
default min_likelihood 0 counts any stated code.

Gate 1: matched mean age >= cohort.gate1_min_mean_age. Design A is the headline
analysis (larger usable n); design B is a sensitivity analysis.

Outputs: out/cohort.csv, out/table1.csv, out/gate1.csv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pulse_ecg import config, io  # noqa: E402

PTBXL_AGE_CAP = 90      # PTB-XL encodes ages > 89 as 300 for privacy


def has_code(codes: dict, wanted, min_lik: float) -> bool:
    return any(c in codes and codes[c] >= min_lik for c in wanted)


def eligible(db: pd.DataFrame, c: dict, min_lik: float) -> pd.DataFrame:
    db = db.copy()
    db["age"] = db["age"].where(db["age"] < 300, PTBXL_AGE_CAP)
    # PTB-XL rhythm statements always carry likelihood 0, so the threshold applies
    # to the diagnostic labels only; a stated sinus code always counts.
    db["sinus"] = db.scp_codes.apply(has_code, args=(c["sinus_codes"], 0))
    db["lae"] = db.scp_codes.apply(has_code, args=(c["primary"]["label_codes"], min_lik))
    db["af"] = db.scp_codes.apply(has_code, args=(c["secondary"]["af_codes"], min_lik))
    db["date"] = pd.to_datetime(db.recording_date)
    root = config.path("ptbxl")
    db["file_present"] = [(root / f"{p}.dat").exists() for p in db.filename_hr]
    return db


def first_per_patient(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(["patient_id", "date"]).groupby("patient_id").head(1)


def design_a(db):
    lae_pat = set(db.loc[db.lae, "patient_id"])
    base = db[db.sinus]
    return (first_per_patient(base[base.lae]),
            first_per_patient(base[~base.patient_id.isin(lae_pat)]))


def design_b(db, timing: str):
    af = db[db.af]
    first_af = af.groupby("patient_id").date.min()
    base = db[db.sinus & ~db.af]
    cases = base[base.patient_id.isin(first_af.index)]
    if timing == "before_first_af":
        cases = cases[cases.date < cases.patient_id.map(first_af)]
    elif timing != "any":
        raise ValueError(f"af_record_timing must be any | before_first_af, got {timing!r}")
    return first_per_patient(cases), first_per_patient(base[~base.patient_id.isin(first_af.index)])


def age_match(cases, controls, tol, k, seed, exact=("sex",)):
    """Greedy k:1 matching: exact on `exact` columns, age within tol (closest
    first), cases in random order."""
    rng = np.random.default_rng(seed)
    pool = controls.assign(_r=rng.random(len(controls)))
    free = pd.Series(True, index=pool.index)
    rows = []
    for ecg_id in rng.permutation(cases.index):
        age = cases.at[ecg_id, "age"]
        ok = free & ((pool.age - age).abs() <= tol)
        for col in exact:
            ok &= pool[col] == cases.at[ecg_id, col]
        cand = pool[ok]
        if cand.empty:
            continue
        pick = cand.assign(_d=(cand.age - age).abs()).sort_values(["_d", "_r"]).index[:k]
        free[pick] = False
        rows.append((ecg_id, "case", ecg_id))
        rows += [(c, "control", ecg_id) for c in pick]
    return pd.DataFrame(rows, columns=["ecg_id", "group", "match_id"])


def build(db, c, design):
    cases, controls = design_a(db) if design == "A" else \
        design_b(db, c["secondary"]["af_record_timing"])
    raw = {"cases_eligible": len(cases), "controls_eligible": len(controls)}
    usable = lambda d: d[d.file_present & (d.age >= c["min_age"])]  # noqa: E731
    cases, controls = usable(cases), usable(controls)
    raw.update(cases_usable=len(cases), controls_usable=len(controls))
    m = age_match(cases, controls, c["age_tol"], c["controls_per_case"], c["seed"],
                  tuple(c["match_exact"]))
    m = m.join(db[["patient_id", "age", "sex", "filename_hr"]], on="ecg_id")
    return m.assign(design=design), raw


def table1(cohort):
    rows = []
    for (d, g), x in cohort.groupby(["design", "group"]):
        rows.append({"design": d, "group": g, "n": len(x),
                     "age_mean": x.age.mean(), "age_sd": x.age.std(),
                     "male_pct": 100 * (x.sex == 0).mean()})      # PTB-XL: 0 = male
    return pd.DataFrame(rows)


if __name__ == "__main__":
    cfg = config.load()
    c = cfg["cohort"]
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-likelihood", type=float, default=c["min_likelihood"])
    args = ap.parse_args()

    db = eligible(io.ptbxl_database(), c, args.min_likelihood)
    print(f"PTB-XL: {len(db)} records, {db.patient_id.nunique()} patients, "
          f"{db.file_present.sum()} with records500 present; "
          f"sinus {db.sinus.sum()}, LAO/LAE {db.lae.sum()}, AFIB {db.af.sum()} "
          f"(min likelihood {args.min_likelihood:g})")

    parts, gate = [], []
    for design in ("A", "B"):
        m, raw = build(db, c, design)
        parts.append(m)
        matched = m[m.group == "case"]
        gate.append({"design": design, "min_likelihood": args.min_likelihood, **raw,
                     "cases_matched": len(matched), "controls_matched": (m.group == "control").sum(),
                     "mean_age_matched": m.age.mean(),
                     "gate1_age_pass": bool(m.age.mean() >= c["gate1_min_mean_age"])})
    cohort = pd.concat(parts, ignore_index=True)
    gate = pd.DataFrame(gate)
    headline = gate.sort_values("cases_matched", ascending=False).design.iloc[0]
    gate["headline"] = gate.design == headline

    out = config.path("out")
    out.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.min_likelihood == c["min_likelihood"] else f"_lik{args.min_likelihood:g}"
    cohort.to_csv(out / f"cohort{suffix}.csv", index=False)
    table1(cohort).to_csv(out / f"table1{suffix}.csv", index=False)
    gate.to_csv(out / f"gate1{suffix}.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20,
                           "display.float_format", "{:.1f}".format):
        print(gate.to_string(index=False))
        print(table1(cohort).to_string(index=False))
