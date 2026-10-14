"""
Record-level QC (pulse_ecg/qc.py) on the Design A cohort, current config.
Produces the cohort-flow numbers in paper/notes/qc_exclusions.md.

    .venv/Scripts/python.exe paper/notes/qc_exclusions.py

QC is decided on the feature row alone; group is merged back only afterwards.
Output: out/qc_designA.csv (one row per record, every record) and
        out/qc_flow_designA.csv (exclusions per rule, all / case / control)
"""

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pulse_ecg import config, io, qc  # noqa: E402


def _one(ecg_id: int) -> dict:
    try:
        row = qc.qc_row(io.read_ptbxl(io.ptbxl_database().loc[ecg_id]), config.load())
    except Exception as exc:  # noqa: BLE001 - becomes a failed delineation, not a drop
        row = {"delineation_reason": f"error: {type(exc).__name__}: {exc}"}
    row["ecg_id"] = ecg_id
    return row


def main():
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort["design"] == "A"]
    with ProcessPoolExecutor() as ex:
        rows = pd.DataFrame(ex.map(_one, a["ecg_id"].tolist(), chunksize=8))
    rules = qc.build_rules()
    res = qc.apply(rows, rules)                                  # blind: no group column
    res = a[["ecg_id", "group", "match_id"]].merge(res, on="ecg_id", how="left")
    assert len(res) == len(a) and res["qc_pass"].notna().all()
    res.to_csv(ROOT / "out" / "qc_designA.csv", index=False)
    f = qc.flow(res, rules, by="group")
    f.to_csv(ROOT / "out" / "qc_flow_designA.csv", index=False)

    wide = f.pivot_table(index="rule", columns="group", values=["any", "first"], sort=False)
    print(wide.to_string())
    kept = res[res["qc_pass"]]
    sets = kept.groupby("match_id")["group"].agg(lambda g: ("case" in set(g), (g == "control").sum()))
    print("\nmatched sets:", res["match_id"].nunique(),
          " case kept:", sum(c for c, _ in sets),
          " case kept + >=1 control:", sum(c and k >= 1 for c, k in sets),
          " case kept + 2 controls:", sum(c and k >= 2 for c, k in sets))


if __name__ == "__main__":
    main()
