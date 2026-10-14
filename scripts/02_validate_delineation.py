"""
02_validate_delineation.py — Gate 2: delineation error against cardiologist
annotations. Paper Table II.

    python scripts/02_validate_delineation.py            # ludb + qtdb
    python scripts/02_validate_delineation.py ludb

LUDB (primary): 12-lead, 500 Hz, 10 s, boundaries annotated per lead. The global
reference for a beat is the earliest onset / latest offset over the leads that
annotate it, which is what the vector-magnitude method measures. LUDB annotates
only the middle of each record, so only reference beats inside the annotated span
are scored; detections outside it are ignored, not counted.

QTDB (secondary): 2-lead, 250 Hz Holter excerpts, one annotation set (q1c) for
a ~30-beat span. That span is cut into 10 s windows, as the pipeline would see a
PTB-XL record. At 250 Hz one sample is 4 ms (LUDB: 2 ms); errors are in ms.

LUDB is split 50/50 by record (gate2.split_seed). Constants in `bounds` are tuned
on the "tune" half only; the Gate 2 figure is the "test" half, sinus rhythm (the
pipeline analyses sinus-rhythm records only). All-rhythm rows are reported too.

Outputs:  out/gate2_beats.csv     one row per scored beat and boundary
          out/gate2_summary.csv   per dataset x boundary: n, mean, SD, MAE, ...
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pulse_ecg import config, io  # noqa: E402
from pulse_ecg.delineate import delineate  # noqa: E402
from pulse_ecg.preprocess import preprocess  # noqa: E402
from pulse_ecg.recording import Recording  # noqa: E402
from pulse_ecg.wavedet import select_multilead  # noqa: E402

BOUNDARIES = ("p_onset", "p_offset", "qrs_onset", "qrs_offset", "t_offset")
MATCH_MS = 150          # reference beat <-> detected R-peak
P_BEFORE_MS = 400       # a P wave belongs to the QRS it precedes by at most this
T_AFTER_MS = 700        # a T wave belongs to the QRS it follows by at most this


def reference_beats(waves: pd.DataFrame, fs: int, rule: str = "extreme",
                    cfg: dict | None = None) -> pd.DataFrame:
    """Per-lead waves -> one global row per beat.

    rule="extreme": earliest onset / latest offset over the annotated leads. This
      is an order statistic over up to 12 independent annotations and inflates
      the span by construction (80 % of LUDB sinus records reach P >= 120 ms).
    rule="knn": the same k-of-N agreement rule the pipeline measures with
      (wavedet.select_multilead, config wavedet.multilead)."""
    ml = (cfg or config.load())["wavedet"]["multilead"]

    def pick(vals, key):
        vals = np.asarray(vals, dtype=float)
        if rule == "extreme" or len(vals) == 1:
            return np.nanmin(vals) if key.endswith("onset") else np.nanmax(vals)
        k = ml["k"][key] if isinstance(ml["k"], dict) else ml["k"]
        return select_multilead(vals, key.endswith("onset"), k, ml["delta_ms"][key] * fs / 1000)
    ms = lambda v: v * fs / 1000                                  # noqa: E731
    qrs = waves[waves.wave == "qrs"].sort_values("peak")
    clusters, cur = [], []
    for _, r in qrs.iterrows():
        if cur and r.peak - cur[0].peak > ms(100):
            clusters.append(cur)
            cur = []
        cur.append(r)
    if cur:
        clusters.append(cur)

    rows = []
    for c in clusters:
        c = pd.DataFrame(c)
        row = {"r_ref": float(c.peak.median()),
               "qrs_onset": pick(c.onset, "qrs_onset"), "qrs_offset": pick(c.offset, "qrs_offset")}
        q_on = row["qrs_onset"] if not np.isnan(row["qrs_onset"]) else row["r_ref"]
        q_off = row["qrs_offset"] if not np.isnan(row["qrs_offset"]) else row["r_ref"]
        p = waves[(waves.wave == "p") & (waves.peak < q_on) & (waves.peak > q_on - ms(P_BEFORE_MS))]
        t = waves[(waves.wave == "t") & (waves.peak > q_off) & (waves.peak < q_off + ms(T_AFTER_MS))]
        row["p_onset"] = pick(p.onset, "p_onset") if len(p) else np.nan
        row["p_offset"] = pick(p.offset, "p_offset") if len(p) else np.nan
        row["t_offset"] = pick(t.offset, "t_offset") if len(t) else np.nan
        row["has_p_ref"] = bool(len(p))
        rows.append(row)
    return pd.DataFrame(rows)


def score(rec: Recording, ref: pd.DataFrame, lo: float, hi: float, offset: int = 0):
    """Delineate `rec`, match to reference beats with r_ref in [lo, hi] (absolute
    samples; `offset` is rec's start in the original signal)."""
    d = delineate(preprocess(rec))
    ref = ref[(ref.r_ref >= lo) & (ref.r_ref <= hi)]
    got = {k: v + offset for k, v in d.beat_points().items()}
    r_abs = d.r_peaks + offset
    tol = MATCH_MS * rec.fs / 1000
    rows = []
    for _, b in ref.iterrows():
        row = {"record": rec.record_id, "r_ref": b.r_ref, "detected": False,
               "has_p_ref": b.has_p_ref, "has_p_got": False, "delin_ok": d.ok}
        if len(r_abs):
            k = int(np.argmin(np.abs(r_abs - b.r_ref)))
            if abs(r_abs[k] - b.r_ref) <= tol:
                row["detected"] = True
                row["has_p_got"] = bool(not np.isnan(got["p_onset"][k]))
                for name in BOUNDARIES:
                    row[name + "_err_ms"] = (got[name][k] - b[name]) * 1000 / rec.fs
                row["p_dur_ref_ms"] = (b.p_offset - b.p_onset) * 1000 / rec.fs
        rows.append(row)
    return rows


def ludb_split(cfg) -> set[str]:
    """Record ids in the tune half. bounds.* may be tuned on these only."""
    g = cfg["gate2"]
    ids = sorted(io.ludb_records(), key=lambda s: int(s.split("/")[-1]))
    perm = np.random.default_rng(g["split_seed"]).permutation(ids)
    return {s.split("/")[-1] for s in perm[:int(len(ids) * g["tune_frac"])]}


def ludb_sinus() -> set[str]:
    """Records whose LUDB rhythm is sinus: the population the pipeline analyses."""
    meta = pd.read_csv(config.path("ludb") / "ludb.csv")
    return {str(i) for i, r in zip(meta.ID, meta.Rhythms) if "sinus" in str(r).lower()}


def run_ludb() -> pd.DataFrame:
    """Scored against both reference constructions (column `reference`)."""
    tune, sinus = ludb_split(config.load()), ludb_sinus()
    rows = []
    for rid in io.ludb_records():
        rec, waves = io.read_ludb(rid)
        lo = np.nanmin(waves[["onset", "peak"]].values)
        hi = np.nanmax(waves[["peak", "offset"]].values)
        for rule in ("extreme", "knn"):
            rows += [dict(r, reference=rule)
                     for r in score(rec, reference_beats(waves, rec.fs, rule), lo, hi)]
    df = pd.DataFrame(rows)
    return df.assign(dataset="ludb", fs=500,
                     half=np.where(df.record.isin(tune), "tune", "test"),
                     sinus=df.record.isin(sinus))


def run_qtdb(window_s: float = 10.0) -> pd.DataFrame:
    rows = []
    for rid in io.qtdb_records():
        rec, waves = io.read_qtdb(rid)
        if waves.empty:
            continue
        ref = reference_beats(waves, rec.fs)
        if ref.empty:                           # no QRS annotated: nothing to score
            print(f"  qtdb {rid}: no annotated QRS, skipped")
            continue
        n = int(window_s * rec.fs)
        start = max(0, int(np.nanmin(waves.peak)) - rec.fs)
        stop = min(rec.n_samples, int(np.nanmax(waves.peak)) + rec.fs)
        for a in range(start, stop, n):
            b = min(a + n, rec.n_samples)
            if b - a < n:                       # last window: slide back to full length
                a = max(0, b - n)
            seg = Recording.from_array(rec.signal[a:b], rec.fs, rec.leads, "qtdb",
                                       rid)
            # score only beats whose whole P-QRS-T fits comfortably inside the window
            rows += score(seg, ref, a + 0.5 * rec.fs, b - 0.7 * rec.fs, offset=a)
    df = pd.DataFrame(rows).drop_duplicates(["record", "r_ref"])
    return df.assign(dataset="qtdb", fs=250, half="external", sinus=np.nan,
                     reference="single annotator")


def summarise(beats: pd.DataFrame, tol: dict, ref_sd: dict) -> pd.DataFrame:
    out = []
    subsets = []
    for (ds, half, ref), g in beats.groupby(["dataset", "half", "reference"]):
        subsets.append((ds, half, ref, "all rhythms", g))
        if g.sinus.notna().any():
            subsets.append((ds, half, ref, "sinus", g[g.sinus.astype(bool)]))
    for ds, half, ref, rhythm, g in subsets:
        det = g.detected.mean()
        pref, pgot = g.has_p_ref, g.has_p_got & g.detected
        for name in BOUNDARIES:
            e = g[name + "_err_ms"].dropna() if name + "_err_ms" in g else pd.Series(dtype=float)
            out.append({
                "dataset": ds, "half": half, "reference": ref, "rhythm": rhythm, "boundary": name, "sample_ms": 1000 / g.fs.iloc[0],
                "n_ref_beats": len(g), "beat_sensitivity": det,
                "p_sensitivity": (pref & pgot).sum() / max(pref.sum(), 1),
                "p_ppv": (pref & pgot).sum() / max(pgot.sum(), 1),
                "n": len(e), "mean_ms": e.mean(), "sd_ms": e.std(),
                # Martinez 2004's statistic: SD within each record, averaged
                "sd_within_record_ms": g.groupby("record")[name + "_err_ms"].std().mean()
                if name + "_err_ms" in g else np.nan,
                "martinez2004_sd_ms": ref_sd.get(name, np.nan),
                "mae_ms": e.abs().mean(), "median_ae_ms": e.abs().median(),
                "p90_ae_ms": e.abs().quantile(0.9),
                "cse_2sd_tol_ms": tol.get(name, np.nan),
                "sd_within_tol": bool(e.std() <= tol.get(name, np.inf)) if len(e) else False,
            })
    return pd.DataFrame(out)


def p_duration_agreement(beats: pd.DataFrame, thr_ms: float) -> pd.DataFrame:
    """P duration is the headline atrial marker, so it is validated as a quantity
    in its own right: per-record median P duration, ours vs reference, and
    agreement on the >= thr_ms call (Cohen's kappa)."""
    rows = []
    for (ds, half, ref), g in beats.groupby(["dataset", "half", "reference"]):
        if ds == "ludb":
            g = g[g.sinus.astype(bool)]
        g = g.dropna(subset=["p_onset_err_ms", "p_offset_err_ms"])
        err = g.p_offset_err_ms - g.p_onset_err_ms
        rec = g.assign(ref=g.p_dur_ref_ms, got=g.p_dur_ref_ms + err).groupby("record")[["ref", "got"]].median()
        e = rec.got - rec.ref
        a, b = rec.ref >= thr_ms, rec.got >= thr_ms
        po = (a == b).mean()
        pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
        rows.append({"dataset": ds, "half": half, "reference": ref,
                     "n_records": len(rec), "n_beats": len(g),
                     "beat_mean_err_ms": err.mean(), "beat_sd_err_ms": err.std(),
                     "record_mean_err_ms": e.mean(), "record_sd_err_ms": e.std(),
                     "record_mae_ms": e.abs().mean(),
                     "ref_prev_ge_thr": a.mean(), "got_prev_ge_thr": b.mean(),
                     "kappa_ge_thr": (po - pe) / (1 - pe) if pe < 1 else np.nan})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    cfg = config.load()
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="*", default=[cfg["gate2"]["primary"], cfg["gate2"]["secondary"]])
    args = ap.parse_args()
    runners = {"ludb": run_ludb, "qtdb": run_qtdb}
    beats = pd.concat([runners[d]() for d in args.datasets], ignore_index=True)
    out = config.path("out")
    out.mkdir(parents=True, exist_ok=True)
    beats.to_csv(out / "gate2_beats.csv", index=False)
    summ = summarise(beats, cfg["gate2"]["cse_2sd_tol_ms"], cfg["gate2"]["martinez2004_sd_ms"])
    summ.to_csv(out / "gate2_summary.csv", index=False)
    pdur = p_duration_agreement(beats, cfg["gate2"]["p_duration_threshold_ms"])
    pdur.to_csv(out / "gate2_p_duration.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20,
                           "display.float_format", "{:.2f}".format):
        print(summ[["dataset", "half", "reference", "rhythm", "boundary", "sample_ms", "n", "beat_sensitivity",
                    "p_sensitivity", "p_ppv", "mean_ms", "sd_ms", "sd_within_record_ms",
                    "martinez2004_sd_ms", "mae_ms", "median_ae_ms", "cse_2sd_tol_ms"]].to_string(index=False))
        print()
        print("P duration (per-record median; kappa on >= threshold):")
        print(pdur.to_string(index=False))
