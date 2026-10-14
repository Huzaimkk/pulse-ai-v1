"""
Feature registry — declare a feature once, derive everything from it.

    @feature(name="p_duration_ms", family="atrial", units="ms", leads=["vm"],
             requires=["p_onset", "p_offset"], normal=(80, 120),
             abnormal=(">=", 120), cite="martinezselles2020bayes",
             desc="Global P-wave duration on the 12-lead vector magnitude")
    def p_duration_ms(ctx):
        return (ctx.p_offset - ctx.p_onset) / ctx.fs * 1000

From that one declaration:
  compute()         the extraction, with a NaN + stated reason whenever a feature
                    cannot run (missing lead, recording too short, failed delineation)
  capability()      which features a lead set / duration can produce (device map)
  table()           Table III (thresholds + citations) and the data dictionary
  sanity()          Gate 3: how much of a feature table lands in the normal range

`abnormal` is an (operator, value) pair rather than a lambda so that it can be
printed in Table III as well as evaluated.
"""

import math
import operator
import re
from dataclasses import dataclass
from typing import Any, Callable, Mapping

import numpy as np
import pandas as pd

from .recording import STANDARD_LEADS, Recording, canonical_lead

FAMILIES = ("atrial", "ventricular", "rate", "quality")
LEAD_TOKENS = {
    "vm": "vector magnitude over the leads present",
    "any": "any single lead",
}
_OPS = {">=": operator.ge, ">": operator.gt, "<=": operator.le,
        "<": operator.lt, "==": operator.eq}
_NAME = re.compile(r"^[a-z][a-z0-9_]*$")


class NotApplicable(Exception):
    """Raised by a feature that is well defined but does not apply to this record
    (e.g. Das narrow-QRS fQRS criteria on a QRS >= 120 ms). The message becomes the
    NaN reason; unlike an error it is never re-raised, even in strict mode."""


@dataclass(frozen=True)
class Threshold:
    op: str
    value: float

    def __post_init__(self):
        if self.op not in _OPS:
            raise ValueError(f"threshold operator must be one of {list(_OPS)}, got {self.op!r}")

    def __call__(self, v) -> bool:
        return bool(_OPS[self.op](v, self.value))

    def __str__(self):
        return f"{self.op} {self.value:g}"


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    func: Callable[["Context"], Any]
    family: str
    units: str
    leads: tuple[str, ...]
    requires: tuple[str, ...]
    min_duration_s: float
    normal: tuple[float, float] | None
    abnormal: Threshold | None
    cite: tuple[str, ...]
    desc: str

    @property
    def reason_key(self) -> str:
        return f"{self.name}_reason"


class Context:
    """What a feature function sees: the recording plus delineation outputs.

    Attribute access falls through to `points` (p_onset, qrs_offset, ...) and
    then to extra keyword values (median beat, vector magnitude, ...)."""

    def __init__(self, rec: Recording, points: Mapping[str, Any] | None = None, **extra):
        self.rec = rec
        self.fs = rec.fs
        self.points = dict(points or {})
        self.extra = extra

    def __getattr__(self, name):
        for store in ("points", "extra"):
            d = self.__dict__.get(store, {})
            if name in d:
                return d[name]
        raise AttributeError(f"context has no {name!r}")

    def has(self, name: str) -> bool:
        v = self.points.get(name, self.extra.get(name))
        if v is None:
            return False
        if np.isscalar(v):
            return not (isinstance(v, float) and math.isnan(v))
        return np.size(v) > 0 and not np.all(np.isnan(np.asarray(v, dtype=float)))


class Registry:
    def __init__(self):
        self._specs: dict[str, FeatureSpec] = {}

    # ------------------------------------------------------------ declaring --
    def feature(self, *, name, family, units, leads, requires=(), min_duration_s=0.0,
                normal=None, abnormal=None, cite=(), desc):
        """Decorator registering a feature function. See module docstring."""
        if not _NAME.match(name):
            raise ValueError(f"feature name {name!r} must be snake_case")
        if name in self._specs:
            raise ValueError(f"feature {name!r} is already registered")
        if family not in FAMILIES:
            raise ValueError(f"{name}: family must be one of {FAMILIES}, got {family!r}")
        leads = tuple(leads)
        if not leads:
            raise ValueError(f"{name}: declare at least one lead (or 'vm' / 'any')")
        leads = tuple(l if l in LEAD_TOKENS else canonical_lead(l) for l in leads)
        bad = [l for l in leads if l not in LEAD_TOKENS and l not in STANDARD_LEADS]
        if bad:
            raise ValueError(f"{name}: unknown lead(s) {bad}")
        if normal is not None:
            lo, hi = normal
            if not lo < hi:
                raise ValueError(f"{name}: normal range must satisfy lo < hi, got {normal}")
            normal = (float(lo), float(hi))
        if abnormal is not None and not isinstance(abnormal, Threshold):
            abnormal = Threshold(*abnormal)
        cite = (cite,) if isinstance(cite, str) else tuple(cite)
        if (normal is not None or abnormal is not None) and not cite:
            raise ValueError(f"{name}: every threshold needs a citation (cite=...)")
        if min_duration_s < 0:
            raise ValueError(f"{name}: min_duration_s must be >= 0")

        def register(func):
            self._specs[name] = FeatureSpec(
                name=name, func=func, family=family, units=units, leads=leads,
                requires=tuple(requires), min_duration_s=float(min_duration_s),
                normal=normal, abnormal=abnormal, cite=cite, desc=desc)
            return func
        return register

    def __iter__(self):
        return iter(self._specs.values())

    def __len__(self):
        return len(self._specs)

    def __getitem__(self, name) -> FeatureSpec:
        return self._specs[name]

    # ------------------------------------------------------------- gating ---
    @staticmethod
    def _lead_reason(spec: FeatureSpec, leads) -> str | None:
        if not leads:
            return "recording has no leads"
        for l in spec.leads:
            if l not in LEAD_TOKENS and l not in leads:
                return f"requires lead {l}, not present"
        return None

    @staticmethod
    def _duration_reason(spec: FeatureSpec, duration_s: float) -> str | None:
        if duration_s < spec.min_duration_s:
            return (f"requires >= {spec.min_duration_s:g} s of signal, "
                    f"recording is {duration_s:g} s")
        return None

    @staticmethod
    def _requires_reason(spec: FeatureSpec, ctx: Context) -> str | None:
        missing = [p for p in spec.requires if not ctx.has(p)]
        if missing:
            return f"requires {', '.join(missing)}, not delineated"
        return None

    def unavailable_reason(self, spec: FeatureSpec, rec: Recording,
                           ctx: Context | None = None) -> str | None:
        ctx = ctx or Context(rec)
        return (self._lead_reason(spec, rec.leads)
                or self._duration_reason(spec, rec.duration_s)
                or self._requires_reason(spec, ctx))

    # ------------------------------------------------------------ compute ---
    def compute(self, rec: Recording, ctx: Context | None = None,
                strict: bool = False) -> dict[str, Any]:
        """One row: {feature: value, feature_reason: "" or why it is NaN}.

        strict=False records an exception as the NaN reason so one bad record
        cannot stop a 21,799-record run; tests pass strict=True so bugs surface."""
        ctx = ctx or Context(rec)
        row: dict[str, Any] = {}
        for spec in self:
            reason = self.unavailable_reason(spec, rec, ctx)
            value: Any = math.nan
            if reason is None:
                try:
                    value = spec.func(ctx)
                    reason = ""
                    if value is None:
                        value, reason = math.nan, "feature returned None"
                except NotApplicable as exc:
                    value, reason = math.nan, str(exc)
                except Exception as exc:          # noqa: BLE001 - reported, not hidden
                    if strict:
                        raise
                    reason = f"error: {type(exc).__name__}: {exc}"
            row[spec.name] = value
            row[spec.reason_key] = reason
        return row

    def is_abnormal(self, name: str, value) -> bool | None:
        spec = self[name]
        if spec.abnormal is None or value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        return spec.abnormal(value)

    # ------------------------------------------------------ derived tables ---
    def capability(self, leads, duration_s: float) -> pd.DataFrame:
        """Which features a lead set and recording length can produce."""
        leads = tuple(canonical_lead(l) for l in leads)
        rows = []
        for spec in self:
            reason = self._lead_reason(spec, leads) or self._duration_reason(spec, duration_s)
            rows.append({"feature": spec.name, "family": spec.family,
                         "available": reason is None, "reason": reason or ""})
        return pd.DataFrame(rows, columns=["feature", "family", "available", "reason"])

    def table(self) -> pd.DataFrame:
        """Data dictionary; rows with a threshold form Table III."""
        rows = [{
            "feature": s.name, "family": s.family, "units": s.units,
            "leads": ", ".join(s.leads), "requires": ", ".join(s.requires),
            "min_duration_s": s.min_duration_s,
            "normal": f"{s.normal[0]:g}–{s.normal[1]:g}" if s.normal else "",
            "abnormal": str(s.abnormal) if s.abnormal else "",
            "cite": ", ".join(s.cite), "description": s.desc,
        } for s in self]
        return pd.DataFrame(rows, columns=["feature", "family", "units", "leads", "requires",
                                           "min_duration_s", "normal", "abnormal", "cite",
                                           "description"])

    def sanity(self, features: pd.DataFrame) -> pd.DataFrame:
        """Gate 3: per feature, the median and the share of values in the normal range."""
        rows = []
        for s in self:
            if s.normal is None or s.name not in features:
                continue
            v = pd.to_numeric(features[s.name], errors="coerce").dropna()
            lo, hi = s.normal
            rows.append({"feature": s.name, "n": len(v),
                         "median": v.median() if len(v) else math.nan,
                         "normal": f"{lo:g}–{hi:g}",
                         "frac_in_normal": ((v >= lo) & (v <= hi)).mean() if len(v) else math.nan})
        return pd.DataFrame(rows, columns=["feature", "n", "median", "normal", "frac_in_normal"])

    def cite_keys(self) -> set[str]:
        return {k for s in self for k in s.cite}


REGISTRY = Registry()
feature = REGISTRY.feature
