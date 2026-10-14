"""
QRS onset vs P offset on the median beat, Design A cohort, current config.
One-off measurement for paper/notes/qrs_onset_overlap.md; not a pipeline stage.
Runs the production context() + REGISTRY.compute (same path as extract()).

    .venv/Scripts/python.exe paper/notes/qrs_onset_overlap.py

Output: out/qrs_onset_overlap_designA.csv (one row per record)
"""

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pulse_ecg import config, io  # noqa: E402
from pulse_ecg.extract import context  # noqa: E402
from pulse_ecg.registry import REGISTRY  # noqa: E402

PTS = ("p_onset", "p_offset", "qrs_onset", "qrs_offset")
FEATS = ("p_axis_deg", "p_duration_ms", "qrs_duration_ms")


def _one(ecg_id: int) -> dict:
    import pulse_ecg.features  # noqa: F401
    out = {"ecg_id": ecg_id}
    try:
        cfg = config.load()
        rec = io.read_ptbxl(io.ptbxl_database().loc[ecg_id])
        ctx = context(rec, cfg)
        fs = ctx.rec.fs
        out["fs"] = fs
        out["method"] = ctx.delineation.method
        out["delineation_reason"] = ctx.delineation.reason
        for k in PTS:
            out[k] = ctx.delineation.points.get(k, np.nan)
        row = REGISTRY.compute(ctx.rec, ctx)
        for f in FEATS:
            out[f] = row.get(f, np.nan)
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def main():
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort["design"] == "A"]
    with ProcessPoolExecutor() as ex:
        rows = list(ex.map(_one, a["ecg_id"].tolist(), chunksize=8))
    res = a[["ecg_id", "group", "match_id"]].merge(pd.DataFrame(rows), on="ecg_id")
    res["gap_ms"] = (res["qrs_onset"] - res["p_offset"]) * 1000 / res["fs"]
    res["pr_ms"] = (res["qrs_onset"] - res["p_onset"]) * 1000 / res["fs"]  # fiducial PR
    res.to_csv(ROOT / "out" / "qrs_onset_overlap_designA.csv", index=False)

    both = res["gap_ms"].notna()
    ov = both & (res["gap_ms"] < 0)
    print(f"n = {len(res)}  errors = {res.get('error', pd.Series(dtype=str)).notna().sum()}"
          f"  both boundaries present = {both.sum()}")
    print(f"overlap overall: {ov.sum()}/{len(res)} = {100 * ov.mean():.2f} % of all records;"
          f" {ov.sum()}/{both.sum()} = {100 * ov[both].mean():.2f} % of records with both")
    for g, sub in res.groupby("group"):
        b, o = sub["gap_ms"].notna(), sub["gap_ms"] < 0
        print(f"  {g:8s} {o.sum()}/{len(sub)} = {100 * o.mean():.2f} %  "
              f"(of {b.sum()} with both: {100 * o[b].mean():.2f} %)")
    gap = res.loc[both, "gap_ms"]
    print(f"gap qrs_onset - p_offset (ms): median {gap.median():.1f}  "
          f"p10 {gap.quantile(.10):.1f}  p90 {gap.quantile(.90):.1f}  "
          f"min {gap.min():.1f}  max {gap.max():.1f}")
    for g, sub in res[both].groupby("group"):
        print(f"  {g:8s} median {sub['gap_ms'].median():.1f}  p90 {sub['gap_ms'].quantile(.9):.1f}")
    print("gap <= 10 ms:", int((gap <= 10).sum()), " gap <= 20 ms:", int((gap <= 20).sum()))
    print("\naffected records:")
    print(res.loc[ov, ["ecg_id", "group", "p_offset", "qrs_onset", "gap_ms", "pr_ms", *FEATS]]
          .to_string(index=False))
    pa = res.loc[both, "p_axis_deg"].dropna()
    print(f"\np_axis_deg all records: median {pa.median():.1f}; "
          f"|axis| > 120: {int((pa.abs() > 120).sum())}/{len(pa)}")


if __name__ == "__main__":
    main()
