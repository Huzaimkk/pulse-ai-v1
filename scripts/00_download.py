"""
00_download.py — fetch PTB-XL (records500 only), LUDB and the QT Database.

Every file is checked against the dataset's own SHA256SUMS.txt, so a finished
run is a verified copy of the published release. Re-running resumes: files
already present with the right hash are skipped.

    python scripts/00_download.py                 # all three
    python scripts/00_download.py ludb qtdb       # just these
"""

import argparse
import hashlib
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pulse_ecg import config  # noqa: E402

CHUNK = 1 << 20


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url: str, dest: Path, expected: str, retries: int = 5) -> int:
    if dest.exists() and sha256(dest) == expected:
        return 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r, open(tmp, "wb") as fh:
                while block := r.read(CHUNK):
                    fh.write(block)
            got = sha256(tmp)
            if got != expected:
                raise IOError(f"sha256 mismatch for {dest.name}")
            tmp.replace(dest)
            return dest.stat().st_size
        except Exception:                                   # noqa: BLE001
            if attempt == retries:
                tmp.unlink(missing_ok=True)
                raise
            time.sleep(2 * attempt)
    return 0


def manifest(base: str) -> list[tuple[str, str]]:
    with urllib.request.urlopen(f"{base}/SHA256SUMS.txt", timeout=60) as r:
        lines = r.read().decode().splitlines()
    return [tuple(reversed(line.split(maxsplit=1))) for line in lines if line.strip()]


def download(name: str, cfg: dict) -> None:
    spec = cfg["download"]["datasets"][name]
    base = f"{cfg['download']['mirror']}/{spec['slug']}"
    root = config.path(name)
    files = [(rel, h) for rel, h in manifest(base)
             if not any(rel.startswith(x) for x in spec["exclude"])]
    print(f"{name}: {len(files)} files -> {root}")

    t0, done, nbytes, failed = time.time(), 0, 0, []
    with ThreadPoolExecutor(max_workers=cfg["download"]["workers"]) as pool:
        futs = {pool.submit(fetch, f"{base}/{rel}", root / rel, h): rel for rel, h in files}
        for fut in as_completed(futs):
            try:
                nbytes += fut.result()
            except Exception as exc:                        # noqa: BLE001
                failed.append(f"{futs[fut]}: {exc}")
            done += 1
            if done % 2000 == 0 or done == len(files):
                el = time.time() - t0
                print(f"  {done:6d}/{len(files)}  {nbytes / 1e6:8.1f} MB  "
                      f"{nbytes / 1e6 / max(el, 1e-9):5.1f} MB/s  {el / 60:5.1f} min", flush=True)
    (root / "SHA256SUMS.txt").write_text(
        "".join(f"{h} {rel}\n" for rel, h in files), encoding="utf-8")
    if failed:
        print(f"  {len(failed)} FAILED (re-run to retry):")
        for f in failed[:20]:
            print("   ", f)
        raise SystemExit(1)
    print(f"  {name}: all {len(files)} files verified\n")


if __name__ == "__main__":
    cfg = config.load()
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="*", default=list(cfg["download"]["datasets"]))
    for name in ap.parse_args().datasets:
        download(name, cfg)
