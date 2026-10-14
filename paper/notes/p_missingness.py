"""
P-wave feature missingness on the Design A cohort at the configured p_span_ms.
One-off measurement for paper/notes/p_missingness.md (Limitations); not a
pipeline stage. Runs the production extract() on every Design A record.

    .venv/Scripts/python.exe paper/notes/p_missingness.py

Output: out/p_missingness_designA.csv (one row per record: value + NaN reason)
"""

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pulse_ecg import config, io  # noqa: E402
from pulse_ecg.extract import extract  # noqa: E402

FEATURES = ("p_duration_ms", "ptfv1_uv_ms", "p_axis_deg")


def _one(ecg_id: int) -> dict:
    db = io.ptbxl_database()
    try:
        row = extract(io.read_ptbxl(db.loc[ecg_id]))
    except Exception as exc:  # noqa: BLE001
        row = {f"{f}_reason": f"error: {type(exc).__name__}: {exc}" for f in FEATURES}
    out = {"ecg_id": ecg_id, "delineation_reason": row.get("delineation_reason", "")}
    for f in FEATURES:
        out[f] = row.get(f, float("nan"))
        out[f"{f}_reason"] = row.get(f"{f}_reason", "")
    return out


def main():
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort["design"] == "A"]
    with ProcessPoolExecutor() as ex:
        rows = list(ex.map(_one, a["ecg_id"].tolist(), chunksize=8))
    res = a[["ecg_id", "group", "match_id"]].merge(pd.DataFrame(rows), on="ecg_id")
    res.insert(0, "p_span_ms", config.load()["wavedet"]["p_span_ms"])
    res.to_csv(ROOT / "out" / "p_missingness_designA.csv", index=False)

    print(f"p_span_ms = {res['p_span_ms'].iat[0]}   n = {len(res)}")
    for f in FEATURES:
        miss = res[f].isna()
        by = miss.groupby(res["group"]).agg(["sum", "count", "mean"])
        print(f"\n{f}: overall {miss.sum()}/{len(res)} = {100 * miss.mean():.1f} %")
        for g, r in by.iterrows():
            print(f"  {g:8s} {int(r['sum'])}/{int(r['count'])} = {100 * r['mean']:.1f} %")
        print("  reasons:")
        print(res.loc[miss].groupby("group")[f"{f}_reason"].value_counts().to_string())


if __name__ == "__main__":
    main()
