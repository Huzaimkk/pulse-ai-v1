"""
POST-HOC REVISION ANALYSIS (not pre-specified in analysis_plan.md).

PTFV1 against LUDB's manual V1 annotations. pulse_ecg/ is not modified; every
number goes through the pipeline's own functions:
  - our V1 P boundaries: wavedet.delineate_lead on the preprocessed V1 lead with
    the pipeline's R-peaks. Its per-record median is asserted to equal
    Delineation.per_lead["V1"], which is the value ptfv1_uv_ms reads.
  - PTFV1 (a) = features.atrial.ptfv1_uv_ms(ctx), the registered feature.
  - PTFV1 (b) = the same median beat and the same terminal_negative(), with the
    window taken from LUDB's V1 P onset/offset (median position relative to R).
Subsets: Gate 2 held-out half, sinus (gate2.split_seed/tune_frac, as in
scripts/02_validate_delineation.py); all sinus records; all 200 records.
Part 4 recomputes PTFV1 depth and duration for the 516 analysed PTB-XL controls.

    .venv/Scripts/python.exe paper/notes/revision_ptfv1_ludb.py

Outputs: out/revision_ptfv1_ludb.csv (summary metrics),
         out/revision_ptfv1_ludb_records.csv, out/revision_ptfv1_ludb_beats.csv,
         out/revision_ptfv1_ptbxl_controls.csv
"""

import importlib.util
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pulse_ecg import config, io  # noqa: E402
from pulse_ecg.extract import context  # noqa: E402
from pulse_ecg.features import atrial  # noqa: E402
from pulse_ecg.registry import NotApplicable  # noqa: E402
from pulse_ecg.wavedet import delineate_lead  # noqa: E402

P_BEFORE_MS = 400            # as Gate 2: a P belongs to the QRS it precedes by <= 400 ms
THR = -4000                  # config atrial.ptfv1_uv_ms abnormal: <= -4000 µV·ms


def _gate2():
    spec = importlib.util.spec_from_file_location(
        "gate2", ROOT / "scripts" / "02_validate_delineation.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _ptf(ctx, on, off):
    """(ptfv1, terminal duration ms, depth µV) on the V1 median beat in [on, off]."""
    try:
        dur, depth = atrial.terminal_negative(atrial._p_window(ctx, "V1", on, off), ctx.fs)
    except NotApplicable:
        return math.nan, math.nan, math.nan
    return -depth * dur, dur, depth


def one_ludb(rid):
    cfg = config.load()
    rec, waves = io.read_ludb(rid)
    rid = str(rid).split("/")[-1]                 # "data/12" -> "12", as Gate 2
    ctx = context(rec, cfg)
    d, fs = ctx.delineation, rec.fs
    ms = 1000.0 / fs
    j = rec.leads.index("V1")
    v1 = delineate_lead(ctx.rec.signal[:, j], fs, d.r_peaks, cfg)
    pl = (d.per_lead or {}).get("V1", {})
    for k in ("p_onset", "p_offset"):              # same code path as the feature
        mine = d.r0 + (float(np.nanmedian(v1[k] - d.r_peaks))
                       if np.any(~np.isnan(v1[k])) else math.nan)
        assert (math.isnan(mine) and math.isnan(pl.get(k, math.nan))) or \
            np.isclose(mine, pl[k]), (rid, k)

    w = waves[waves.lead == "V1"]
    lo, hi = np.nanmin(w[["onset", "peak"]].values), np.nanmax(w[["peak", "offset"]].values)
    lp = w[w.wave == "p"]
    beats = []
    for b, r in enumerate(d.r_peaks):
        if not lo <= r <= hi:
            continue
        cand = lp[(lp.peak < r) & (lp.peak > r - P_BEFORE_MS / ms)]
        ref = cand.iloc[-1] if len(cand) else None
        ours = not (np.isnan(v1["p_onset"][b]) or np.isnan(v1["p_offset"][b]))
        row = {"record": rid, "beat": b, "r": int(r), "ludb_has_p": ref is not None,
               "ours_has_p": ours}
        if ref is not None:
            row["ludb_on_rel"] = (ref.onset - r) * ms
            row["ludb_off_rel"] = (ref.offset - r) * ms
            if ours:
                row["err_on_ms"] = (v1["p_onset"][b] - ref.onset) * ms
                row["err_off_ms"] = (v1["p_offset"][b] - ref.offset) * ms
        beats.append(row)
    bt = pd.DataFrame(beats)

    try:
        a = atrial.ptfv1_uv_ms(ctx)
    except NotApplicable:
        a = math.nan
    on_o, off_o = pl.get("p_onset", math.nan), pl.get("p_offset", math.nan)
    a2, dur_a, dep_a = _ptf(ctx, on_o, off_o)
    assert (math.isnan(a) and math.isnan(a2)) or np.isclose(a, a2), rid
    if len(bt) and bt.ludb_has_p.any():
        on_l = d.r0 + bt.ludb_on_rel.median() / ms
        off_l = d.r0 + bt.ludb_off_rel.median() / ms
    else:
        on_l = off_l = math.nan
    b_, dur_b, dep_b = _ptf(ctx, on_l, off_l)
    rec_row = {"record": rid, "n_beats_in_span": len(bt),
               "ludb_has_p_v1": bool(len(lp)), "ours_has_p_v1": not math.isnan(on_o) and not math.isnan(off_o),
               "ours_on_rel_ms": (on_o - d.r0) * ms, "ours_off_rel_ms": (off_o - d.r0) * ms,
               "ludb_on_rel_ms": (on_l - d.r0) * ms, "ludb_off_rel_ms": (off_l - d.r0) * ms,
               "ptfv1_ours": a, "term_dur_ours_ms": dur_a, "depth_ours_uv": dep_a,
               "ptfv1_ludb": b_, "term_dur_ludb_ms": dur_b, "depth_ludb_uv": dep_b}
    rec_row["err_on_ms"] = rec_row["ours_on_rel_ms"] - rec_row["ludb_on_rel_ms"]
    rec_row["err_off_ms"] = rec_row["ours_off_rel_ms"] - rec_row["ludb_off_rel_ms"]
    return rec_row, bt


def one_ptbxl(ecg_id):
    cfg = config.load()
    ctx = context(io.read_ptbxl(io.ptbxl_database().loc[ecg_id]), cfg)
    pl = (ctx.delineation.per_lead or {}).get("V1", {})
    on, off = pl.get("p_onset", math.nan), pl.get("p_offset", math.nan)
    p, dur, dep = _ptf(ctx, on, off)
    return {"ecg_id": ecg_id, "ptfv1_recomputed": p, "term_dur_ms": dur, "depth_uv": dep,
            "v1_p_dur_ms": (off - on) * 1000 / ctx.fs}


# ------------------------------------------------------------- statistics ---
def ba(x):
    x = pd.Series(x).dropna()
    m, s = x.mean(), x.std()
    return {"n": len(x), "bias": m, "sd": s, "loa_lo": m - 1.96 * s, "loa_hi": m + 1.96 * s,
            "mae": x.abs().mean()}


def icc21(x, y):
    """ICC(2,1): two-way random effects, absolute agreement, single measurement."""
    Y = np.column_stack([x, y]).astype(float)
    n, k = Y.shape
    gm = Y.mean()
    msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
    msc = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    sse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2).sum()
    mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)


def kappa(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe) if pe < 1 else math.nan


def summarise(rec, bt, subset):
    out = []

    def add(section, metric, level, **v):
        out.append({"analysis": "POST-HOC revision", "subset": subset, "section": section,
                    "metric": metric, "level": level, **v})

    for k in ("on", "off"):
        add("1 V1 boundary", f"P {k}set error (ours - LUDB), ms", "per-record",
            **ba(rec[f"err_{k}_ms"]))
        add("1 V1 boundary", f"P {k}set error (ours - LUDB), ms", "per-beat",
            **ba(bt[f"err_{k}_ms"]) if "err_on_ms" in bt else {})
    for lev, df, lc, oc in (("per-record", rec, "ludb_has_p_v1", "ours_has_p_v1"),
                            ("per-beat", bt, "ludb_has_p", "ours_has_p")):
        for lv in (True, False):
            s = df[df[lc] == lv]
            add("1 V1 detection", f"ours detects V1 P | LUDB {'has' if lv else 'has no'} V1 P",
                lev, n=len(s), k=int(s[oc].sum()), pct=100 * s[oc].mean() if len(s) else math.nan)

    both = rec.dropna(subset=["ptfv1_ours", "ptfv1_ludb"])
    if len(both) > 2:
        d = both.ptfv1_ours - both.ptfv1_ludb
        add("2 PTFV1 agreement", "PTFV1 ours - LUDB-boundaries, µV·ms", "per-record", **ba(d))
        add("2 PTFV1 agreement", "Spearman r", "per-record", n=len(both),
            value=stats.spearmanr(both.ptfv1_ours, both.ptfv1_ludb)[0])
        add("2 PTFV1 agreement", "ICC(2,1)", "per-record", n=len(both),
            value=icc21(both.ptfv1_ours, both.ptfv1_ludb))
        ao, al = both.ptfv1_ours <= THR, both.ptfv1_ludb <= THR
        add("2 PTFV1 agreement", f"Cohen kappa, abnormal <= {THR}", "per-record", n=len(both),
            value=kappa(ao, al), k=int((ao & al).sum()),
            note=f"ours abn {ao.sum()}, LUDB abn {al.sum()}, both {(ao & al).sum()}, "
                 f"ours only {(ao & ~al).sum()}, LUDB only {(~ao & al).sum()}")
        for c, lab in (("term_dur", "terminal negative duration, ms"),
                       ("depth", "terminal negative depth, µV")):
            u = "_ms" if c == "term_dur" else "_uv"
            add("2 PTFV1 components", f"{lab} ours - LUDB", "per-record",
                **ba(both[f"{c}_ours{u}"] - both[f"{c}_ludb{u}"]))
    for col, lab in (("ptfv1_ours", "ours"), ("ptfv1_ludb", "LUDB boundaries")):
        v = rec[col].dropna()
        add("3 % abnormal", f"PTFV1 <= {THR}, {lab}", "per-record", n=len(v),
            k=int((v <= THR).sum()), pct=100 * (v <= THR).mean(),
            median=v.median(), q1=v.quantile(.25), q3=v.quantile(.75))
    return out


def main():
    cfg = config.load()
    g2 = _gate2()
    tune, sinus = g2.ludb_split(cfg), g2.ludb_sinus()
    ids = io.ludb_records()
    with ProcessPoolExecutor() as ex:
        res = list(ex.map(one_ludb, ids, chunksize=4))
    rec = pd.DataFrame([r for r, _ in res])
    bt = pd.concat([b for _, b in res if len(b)], ignore_index=True)
    rec["half"] = np.where(rec.record.isin(tune), "tune", "test")
    rec["sinus"] = rec.record.isin(sinus)
    bt = bt.merge(rec[["record", "half", "sinus"]], on="record")
    test_ids = sorted(rec[(rec.half == "test") & rec.sinus].record, key=int)

    rows = []
    for name, m in (("Gate2 held-out half, sinus", (rec.half == "test") & rec.sinus),
                    ("all sinus", rec.sinus), ("all 200 records", rec.sinus | ~rec.sinus)):
        ids_ = set(rec[m].record)
        rows += summarise(rec[m], bt[bt.record.isin(ids_)], name)

    # 3/4. PTB-XL analysed controls: stored values + recomputed components
    qc = pd.read_csv(ROOT / "out" / "qc_designA.csv")
    ctl = qc[qc.qc_pass & (qc.group == "control")][["ecg_id", "ptfv1_uv_ms"]]
    with ProcessPoolExecutor() as ex:
        comp = pd.DataFrame(list(ex.map(one_ptbxl, ctl.ecg_id.tolist(), chunksize=8)))
    ctl = ctl.merge(comp, on="ecg_id")
    same = np.isclose(ctl.ptfv1_uv_ms, ctl.ptfv1_recomputed, equal_nan=True)
    v = ctl.ptfv1_uv_ms.dropna()
    calls = pd.read_csv(ROOT / "out" / "group_threshold_calls_designA.csv")
    stored = calls[calls.feature == "ptfv1_uv_ms"].iloc[0]
    for metric, val in (
            ("n with PTFV1", len(v)),
            ("stored % <= -4000 (group_threshold_calls_designA.csv)", stored.control_pct),
            ("recomputed % <= -4000", 100 * (v <= THR).mean()),
            ("% < -4000 (strict)", 100 * (v < THR).mean()),
            ("n exactly -4000", int((v == THR).sum())),
            ("n > 0 (impossible by definition)", int((v > 0).sum())),
            ("n == 0 (no terminal negative component)", int((v == 0).sum())),
            ("PTFV1 median, µV·ms", v.median()), ("PTFV1 Q1", v.quantile(.25)),
            ("PTFV1 Q3", v.quantile(.75)),
            ("recomputed == stored (records)", int(same.sum())),
            ("terminal duration median, ms", ctl.term_dur_ms.median()),
            ("terminal duration Q1", ctl.term_dur_ms.quantile(.25)),
            ("terminal duration Q3", ctl.term_dur_ms.quantile(.75)),
            ("terminal depth median, µV", ctl.depth_uv.median()),
            ("terminal depth Q1", ctl.depth_uv.quantile(.25)),
            ("terminal depth Q3", ctl.depth_uv.quantile(.75)),
            ("V1 single-lead P duration median, ms", ctl.v1_p_dur_ms.median()),
            ("V1 single-lead P duration Q1", ctl.v1_p_dur_ms.quantile(.25)),
            ("V1 single-lead P duration Q3", ctl.v1_p_dur_ms.quantile(.75))):
        rows.append({"analysis": "POST-HOC revision", "subset": "PTB-XL analysed controls",
                     "section": "4 control PTFV1 check", "metric": metric, "level": "per-record",
                     "value": val})
    # LUDB components for the same comparison (sinus, our boundaries)
    s = rec[rec.sinus]
    for metric, col in (("terminal duration median, ms", "term_dur_ours_ms"),
                        ("terminal depth median, µV", "depth_ours_uv")):
        rows.append({"analysis": "POST-HOC revision", "subset": "all sinus",
                     "section": "4 control PTFV1 check", "metric": metric + " (LUDB, ours)",
                     "level": "per-record", "value": s[col].median()})

    summ = pd.DataFrame(rows)
    out = ROOT / "out"
    summ.to_csv(out / "revision_ptfv1_ludb.csv", index=False)
    rec.to_csv(out / "revision_ptfv1_ludb_records.csv", index=False)
    bt.to_csv(out / "revision_ptfv1_ludb_beats.csv", index=False)
    ctl.to_csv(out / "revision_ptfv1_ptbxl_controls.csv", index=False)

    print("POST-HOC revision: PTFV1 vs LUDB V1 annotations\n")
    print(f"Gate 2 held-out sinus records (n={len(test_ids)}):", ", ".join(test_ids), "\n")
    cols = ["subset", "section", "metric", "level", "n", "k", "pct", "bias", "sd", "loa_lo",
            "loa_hi", "mae", "value", "median", "q1", "q3", "note"]
    with pd.option_context("display.width", 320, "display.max_columns", 30,
                           "display.max_colwidth", 80, "display.max_rows", 500):
        print(summ.reindex(columns=cols).to_string(index=False,
              float_format=lambda x: f"{x:.4g}", na_rep=""))


if __name__ == "__main__":
    main()
