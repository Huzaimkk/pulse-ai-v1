"""
PULSE AI — ECG Track: Full PTB-XL parallel extraction runner
=============================================================
Runs the biomarker extractor over all 21,799 PTB-XL records.

WHY THIS FILE EXISTS
--------------------
The bottleneck in this pipeline is NOT matrix multiplication. It is:
    (a) reading 21,799 WFDB files from disk, and
    (b) running scipy filtering + delineation on each one.
Both are CPU- and I/O-bound, and every record is independent of every other.

That makes this an *embarrassingly parallel* problem, and the correct fix is
process-level parallelism across records. Concretely:

  - A GPU does not help. There is no large tensor operation here.
  - A TPU does not help either, for the same reason. TPUs accelerate the dense
    linear algebra inside neural network training. scipy.signal.filtfilt,
    wfdb.rdrecord and tree-based models do not touch a TPU at all; it would sit
    idle while one CPU core does all the work.
  - Quantum computing does not help. There is no quantum algorithm for
    Butterworth filtering, R-peak detection, or gradient boosting on tabular
    data, and running a quantum simulator on classical hardware is slower than
    the classical code it simulates, not faster.

What actually delivers the speed-up, in order of effect:
  1. multiprocessing across records        (~N_cores x)
  2. Parquet instead of CSV for the output (~5-10x smaller, much faster reads)
  3. checkpointing, so a crash at record 19,000 does not cost the whole run

Usage
-----
    python 04_run_full_ptbxl.py --root /content/drive/MyDrive/ptb-xl \\
                                --out out/biomarkers_full.parquet
    python 04_run_full_ptbxl.py --benchmark        # measure without PTB-XL
    python 04_run_full_ptbxl.py --root ... --limit 200   # smoke test first
"""

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

_ext = import_module("02_extract_biomarkers") if Path(
    Path(__file__).parent / "02_extract_biomarkers.py").exists() else None


# ---------------------------------------------------------------------------
# WORKER
# ---------------------------------------------------------------------------
def _process_one(args):
    """Runs in a separate process. Must be importable at module level."""
    ecg_id, rel_path, root, meta = args
    try:
        import wfdb
        rec = wfdb.rdrecord(str(Path(root) / rel_path))
        bm = _ext.extract(rec.p_signal, record_name=str(ecg_id), fs=int(rec.fs))
        row = asdict(bm)
    except Exception as exc:                                  # noqa: BLE001
        row = asdict(_ext.Biomarkers(record=str(ecg_id),
                                     qc_reason=f"error: {type(exc).__name__}: {exc}"))
    row.update(meta)
    return row


# ---------------------------------------------------------------------------
# RUNNER
# ---------------------------------------------------------------------------
def run(root, out_path, workers=None, limit=None, chunk=500, use_hr=True,
        cohort_csv=None):
    root = Path(root)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    workers = workers or max(1, (os.cpu_count() or 2))

    # ---- which records -----------------------------------------------------
    if cohort_csv:
        db = pd.read_csv(cohort_csv, index_col=0)
    else:
        import ast
        db = pd.read_csv(root / "ptbxl_database.csv", index_col="ecg_id")
        db["scp_codes"] = db["scp_codes"].apply(ast.literal_eval)
    col = "filename_hr" if use_hr else "filename_lr"
    if limit:
        db = db.iloc[:limit]

    # ---- resume from checkpoint -------------------------------------------
    done = set()
    ckpt = out_path.with_suffix(".partial.parquet")
    if ckpt.exists():
        prev = pd.read_parquet(ckpt)
        done = set(prev["record"].astype(str))
        print(f"Resuming: {len(done)} records already done in {ckpt.name}")
    todo = [i for i in db.index if str(i) not in done]

    jobs = [
        (i, db.at[i, col], str(root),
         {"age": db.at[i, "age"] if "age" in db.columns else np.nan,
          "sex": db.at[i, "sex"] if "sex" in db.columns else np.nan,
          "patient_id": db.at[i, "patient_id"] if "patient_id" in db.columns else np.nan,
          "strat_fold": db.at[i, "strat_fold"] if "strat_fold" in db.columns else np.nan,
          "group": db.at[i, "group"] if "group" in db.columns else ""})
        for i in todo
    ]

    print(f"{len(jobs)} records to process on {workers} worker processes\n")
    rows, t0, last = [], time.time(), time.time()

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_process_one, j): j[0] for j in jobs}
        for n, fut in enumerate(as_completed(futures), 1):
            rows.append(fut.result())

            if n % chunk == 0 or n == len(jobs):
                elapsed = time.time() - t0
                rate = n / elapsed
                eta = (len(jobs) - n) / rate if rate else 0
                print(f"  {n:6d}/{len(jobs)}  {rate:6.1f} rec/s  "
                      f"elapsed {elapsed/60:5.1f} min  ETA {eta/60:5.1f} min")
                # checkpoint
                df = pd.DataFrame(rows)
                if ckpt.exists():
                    df = pd.concat([pd.read_parquet(ckpt), df], ignore_index=True)
                df.to_parquet(ckpt, index=False)
                rows = []
                last = time.time()

    final = pd.read_parquet(ckpt) if ckpt.exists() else pd.DataFrame(rows)
    if rows:
        final = pd.concat([final, pd.DataFrame(rows)], ignore_index=True)
    final.to_parquet(out_path, index=False)
    ckpt.unlink(missing_ok=True)

    print(f"\nDone in {(time.time()-t0)/60:.1f} min -> {out_path}")
    _summary(final)
    return final


def _summary(df):
    n = len(df)
    print(f"\n{'':-<52}")
    print(f"  records           {n}")
    if "qc_pass" in df:
        print(f"  QC pass           {int(df.qc_pass.sum())} ({df.qc_pass.mean()*100:.1f}%)")
    if "p_detected" in df:
        print(f"  P wave detected   {int(df.p_detected.sum())} ({df.p_detected.mean()*100:.1f}%)")
    if "patient_id" in df:
        print(f"  unique patients   {df.patient_id.nunique()}")
    if "iab_class" in df:
        vc = df.iab_class.value_counts().to_dict()
        print(f"  IAB class         {vc}")
    if "qc_reason" in df:
        bad = df.loc[df.qc_reason.astype(str).str.len() > 0, "qc_reason"]
        if len(bad):
            print("\n  top QC failure reasons:")
            for r, c in bad.value_counts().head(6).items():
                print(f"    {c:5d}  {str(r)[:70]}")
    print(f"{'':-<52}")


# ---------------------------------------------------------------------------
# BENCHMARK — measure real throughput without needing PTB-XL
# ---------------------------------------------------------------------------
def _bench_workload(seed):
    """
    Stand-in with the same computational profile as extract(): 12-channel
    zero-phase notch + bandpass, R-peak detection, median beat, and boundary
    search on the vector magnitude. Used only when 02_extract_biomarkers.py is
    not importable, so the benchmark still reports honest per-record cost.
    """
    from scipy import signal as sps
    rng = np.random.default_rng(seed)
    fs, n = 500, 5000
    t = np.arange(n) / fs
    hr = 60 + (seed % 40)
    sig = np.zeros((n, 12))
    for k in range(int(10 * hr / 60)):
        c = 0.5 + k * 60.0 / hr
        if c > 9.4:
            break
        for li in range(12):
            s = (0.4 + 0.1 * li)
            sig[:, li] += s * np.exp(-0.5 * ((t - c) / 0.014) ** 2)
            sig[:, li] += 0.12 * s * np.exp(-0.5 * ((t - c + 0.16) / 0.027) ** 2)
            sig[:, li] += 0.30 * s * np.exp(-0.5 * ((t - c - 0.28) / 0.055) ** 2)
    sig += rng.normal(0, 0.005, sig.shape)

    b, a = sps.iirnotch(50 / (fs / 2), 30)
    sig = sps.filtfilt(b, a, sig, axis=0)
    b, a = sps.butter(4, [0.5 / (fs / 2), 40 / (fs / 2)], btype="band")
    sig = sps.filtfilt(b, a, sig, axis=0)

    lead2 = sig[:, 1]
    bb, aa = sps.butter(3, [5 / (fs / 2), 15 / (fs / 2)], btype="band")
    x = sps.filtfilt(bb, aa, lead2)
    x = np.diff(x, prepend=x[0]) ** 2
    x = np.convolve(x, np.ones(75) / 75, mode="same")
    rp, _ = sps.find_peaks(x, distance=125, height=np.percentile(x, 90))

    stack = [sig[r - 200:r + 250, :] for r in rp if r >= 200 and r + 250 <= n]
    if len(stack) >= 3:
        beat = np.median(np.stack(stack), axis=0)
        vm = np.sqrt(np.mean(beat ** 2, axis=1))
        base = np.median(vm)
        pk = int(np.argmax(vm))
        thr = base + 0.05 * (vm[pk] - base)
        lo = pk
        while lo > 0 and vm[lo] > thr:
            lo -= 1
        sps.find_peaks(vm[:pk], prominence=1e-6)
    return seed


def _bench_one(seed):
    if _ext is not None:
        sig = _ext._synthetic_ecg(hr=60 + (seed % 40))
        return asdict(_ext.extract(sig, record_name=str(seed)))
    return _bench_workload(seed)


def benchmark(n=24, workers=None):
    workers = workers or max(1, (os.cpu_count() or 2))
    src = "full extractor" if _ext is not None else "equivalent-cost workload"
    print(f"Benchmark: {n} synthetic 12-lead records @500 Hz, "
          f"{os.cpu_count()} cores ({src})\n")

    # Warm up first: SciPy caches filter designs and the first call in a
    # process pays import and allocation costs that would otherwise be charged
    # entirely to the serial run and fake a superlinear speed-up.
    for i in range(3):
        _bench_one(i)

    t0 = time.time()
    for i in range(n):
        _bench_one(i)
    serial = time.time() - t0
    print(f"  serial      {serial:6.2f} s   {n/serial:6.2f} rec/s")

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        list(pool.map(_bench_one, range(n)))
    par = time.time() - t0
    print(f"  {workers} workers   {par:6.2f} s   {n/par:6.2f} rec/s   "
          f"speed-up {serial/par:.2f}x")

    per_rec = par / n
    print(f"\n  per record  {per_rec*1000:.0f} ms")
    for cores, label in ((2, "this sandbox"), (4, "typical laptop"),
                         (8, "Colab CPU / good laptop"), (16, "workstation")):
        scaled = per_rec * (workers / cores)
        total = 21799 * scaled / 60
        print(f"  21,799 records on {cores:2d} cores ({label:23s}) "
              f"~{total:6.1f} min")
    print("\n  NOTE: real PTB-XL records add disk I/O per file, so expect "
          "roughly\n  1.3-2x these figures on first run (OS page cache makes "
          "re-runs faster).")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=str, help="PTB-XL root directory")
    ap.add_argument("--out", type=str, default="out/biomarkers_full.parquet")
    ap.add_argument("--cohort", type=str, default=None,
                    help="cohort.csv from 01_build_cohort.py (subset instead of all)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--lr", action="store_true",
                    help="use 100 Hz records -- NOT recommended, see roadmap 2.3")
    ap.add_argument("--benchmark", action="store_true")
    a = ap.parse_args()

    if a.benchmark or not a.root:
        benchmark(workers=a.workers)
    else:
        run(a.root, a.out, workers=a.workers, limit=a.limit,
            use_hr=not a.lr, cohort_csv=a.cohort)
