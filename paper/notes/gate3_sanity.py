"""
Gate 3 sanity check and post-QC cohort characteristics (Table I draft), Design A.
Reads existing outputs only; runs no extraction, model or outcome comparison.

    .venv/Scripts/python.exe paper/notes/gate3_sanity.py

Inputs:  out/qc_designA.csv (qc_exclusions.py), out/cohort.csv
Outputs: out/gate3_sanity_designA.csv, out/gate3_offsets_designA.csv,
         out/table1_postqc_designA.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pulse_ecg.features  # noqa: E402,F401  (registers every feature)
from pulse_ecg import config  # noqa: E402
from pulse_ecg.registry import REGISTRY  # noqa: E402

# CLAUDE_1.md Gate 3, plus mean HR as requested
PUBLISHED = {"pr_ms": (120, 200), "qrs_duration_ms": (70, 110), "qt_ms": (350, 450),
             "p_duration_ms": (80, 120), "tpe_qt_ratio": (0.15, 0.25),
             "mean_hr_bpm": (60, 100)}

# Method- and age-matched comparators for the multilead (global) durations.
# The single-lead caliper ranges above understate simultaneous-lead durations
# (kligfield2007aha1). Medians from rijnbeek2014normal, lead-independent, age 60-89.
MATCHED = {"p_duration_ms": {"pooled": (115, 125), "men": (120, 121), "women": (114, 118)},
           "qrs_duration_ms": {"pooled": (95, 100), "men": (98, 101), "women": (92, 92)}}


def summarise(name, v, declared, published):
    v = pd.to_numeric(v, errors="coerce").dropna()
    rng = published or declared
    q1, med, q3 = (v.quantile([.25, .5, .75]) if len(v) else [np.nan] * 3)
    row = {"feature": name, "n": len(v), "median": med, "q1": q1, "q3": q3,
           "declared": f"{declared[0]:g}-{declared[1]:g}" if declared else "",
           "published": f"{published[0]:g}-{published[1]:g}" if published else ""}
    if rng and len(v):
        lo, hi = rng
        row |= {"pct_below": 100 * (v < lo).mean(), "pct_above": 100 * (v > hi).mean(),
                "pct_outside": 100 * ((v < lo) | (v > hi)).mean(),
                "verdict": "PASS" if lo <= med <= hi else "FAIL"}
    else:
        row["verdict"] = "no range"
    return row


def smd_cont(a, b):
    return (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)


def smd_bin(a, b):
    p, q = a.mean(), b.mean()
    return (p - q) / np.sqrt((p * (1 - p) + q * (1 - q)) / 2)


def table1(df):
    ca, co = df[df.group == "case"], df[df.group == "control"]
    rows = []
    for label, x in (("case", ca), ("control", co)):
        rows.append({"group": label, "n": len(x), "age_mean": x.age.mean(),
                     "age_sd": x.age.std(ddof=1), "male_n": int((x.sex == 0).sum()),
                     "male_pct": 100 * (x.sex == 0).mean()})
    rows.append({"group": "SMD", "age_mean": smd_cont(ca.age, co.age),
                 "male_pct": smd_bin((ca.sex == 0).astype(float), (co.sex == 0).astype(float))})
    return pd.DataFrame(rows)


def main():
    fs = config.load()["ptbxl"]["fs"]
    qc = pd.read_csv(ROOT / "out" / "qc_designA.csv")
    kept = qc[qc["qc_pass"]].copy()
    assert len(kept) == 775, len(kept)

    # PR is not a registered feature. Derived from the same median-beat points that
    # p_duration_ms uses: PR = (qrs_onset - p_onset)/fs = P duration + PR segment.
    kept["pr_ms"] = kept["p_duration_ms"] + (kept["qrs_onset"] - kept["p_offset"]) * 1000 / fs

    rows = []
    for s in REGISTRY:
        if s.name == "rr_intervals_ms":          # list-valued, summarised by mean_hr_bpm
            continue
        rows.append(summarise(s.name, kept[s.name], s.normal, PUBLISHED.get(s.name)))
    rows.append(summarise("pr_ms (derived, not registered)", kept["pr_ms"], None,
                          PUBLISHED["pr_ms"]))
    g3 = pd.DataFrame(rows)
    g3.to_csv(ROOT / "out" / "gate3_sanity_designA.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 99):
        print(g3.round(3).to_string())      # pooled only: no outcome split by group

    # Residual offset vs the matched comparators: pooled, and by sex at age 60-89.
    # Pooled over case/control; stratified by sex only, never by group.
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort.design == "A"]
    demo = kept.merge(a[["ecg_id", "age", "sex"]], on="ecg_id")
    assert len(demo) == 775
    strata = {"pooled": demo,
              "men": demo[(demo.sex == 0) & demo.age.between(60, 89)],
              "women": demo[(demo.sex == 1) & demo.age.between(60, 89)]}
    off = []
    for feat, refs in MATCHED.items():
        for stratum, (lo, hi) in refs.items():
            v = strata[stratum][feat].dropna()
            off.append({"feature": feat, "stratum": stratum, "n": len(v),
                        "median": v.median(), "ref_lo": lo, "ref_hi": hi,
                        "offset_lo": v.median() - hi, "offset_hi": v.median() - lo})
    off = pd.DataFrame(off)
    off.to_csv(ROOT / "out" / "gate3_offsets_designA.csv", index=False)
    print(off.to_string())

    post = a[a.ecg_id.isin(kept.ecg_id)]
    assert (post.group == "case").sum() == 259 and (post.group == "control").sum() == 516
    t_pre, t_post = table1(a).assign(stage="pre-QC"), table1(post).assign(stage="post-QC")
    t = pd.concat([t_pre, t_post])
    t.to_csv(ROOT / "out" / "table1_postqc_designA.csv", index=False)
    print(t.round(4).to_string())


if __name__ == "__main__":
    main()
