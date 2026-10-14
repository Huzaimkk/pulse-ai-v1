"""
The one entry point: Recording -> one flat row of features.

    row = extract(rec)          # PTB-XL, LUDB, QTDB or device: same call

Pipeline: preprocess (notch + band-pass) -> delineate (wavedet + multilead rule,
validated at Gate 2) -> Context -> REGISTRY.compute. Features see
the filtered recording as ctx.rec and the unfiltered one as ctx.raw (quality
indices must see the drift and mains that filtering removes).
"""

from typing import Any

from . import config
from .delineate import delineate
from .preprocess import preprocess
from .recording import Recording
from .registry import REGISTRY, Context

META_KEYS = ("patient_id", "age", "sex")


def context(rec: Recording, cfg: dict | None = None) -> Context:
    """Preprocess and delineate `rec`; everything a feature function may read."""
    cfg = cfg or config.load()
    filt = preprocess(rec, cfg)
    d = delineate(filt, cfg)
    return Context(filt, d.points, raw=rec, delineation=d, cfg=cfg, **d.context_extra())


def extract(rec: Recording, cfg: dict | None = None, strict: bool = False) -> dict[str, Any]:
    """One row: identity, delineation status, {feature, feature_reason} for every
    registered feature. strict=True re-raises feature bugs."""
    import pulse_ecg.features  # noqa: F401  (registers every feature)

    cfg = cfg or config.load()
    ctx = context(rec, cfg)
    d = ctx.delineation
    row: dict[str, Any] = {
        "record_id": rec.record_id, "source": rec.source, "fs": rec.fs,
        "n_leads": rec.n_leads, "duration_s": rec.duration_s,
        "delineation_method": d.method, "delineation_reason": d.reason,
        **{k: rec.meta[k] for k in META_KEYS if k in rec.meta},
    }
    row.update(REGISTRY.compute(ctx.rec, ctx, strict=strict))
    return row
