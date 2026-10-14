"""
Record-level quality control: pass/fail plus every rule a record failed.

    rules = build_rules()                 # from REGISTRY + config `qc:`
    res = evaluate(row, rules)            # row = qc_row(rec) or any extract()-style row
    res.passed, res.failed                # [("range:qrs_ms", "qrs_duration_ms 231 ms ..."), ...]
    table = apply(rows_df, rules)         # one output row per input row, never fewer

Nothing is dropped silently: `apply` returns every input record with qc_pass and
the names and details of every failed rule, and `flow` counts exclusions per rule.

evaluate() sees only the feature row. It never takes case/control status, so the
rules are applied blind to group.

Rules, in the order `flow` attributes a record to its first failure:
  delineation            delineation failed (extract()'s delineation_reason)
  missing:<name>         a required feature or fiducial is NaN; its reason is carried
  qrs_onset_before_p_offset
  range:<qc key>         a feature outside the config `qc:` plausibility range
"""

import math
from dataclasses import dataclass, field
from typing import Any, Mapping

import pandas as pd

from . import config
from .registry import REGISTRY, Registry

# config `qc:` key -> registered feature it bounds
QC_FEATURES = {"hr": "mean_hr_bpm", "qrs_ms": "qrs_duration_ms",
               "qt_ms": "qt_ms", "p_ms": "p_duration_ms"}
# median-beat fiducials the ordering rule reads (delineation.points, samples)
FIDUCIALS = ("p_offset", "qrs_onset")
OVERLAP_RULE = "qrs_onset_before_p_offset"


@dataclass(frozen=True)
class RangeRule:
    key: str            # config qc key, e.g. "qrs_ms"
    feature: str        # registered feature name
    lo: float
    hi: float
    units: str

    @property
    def name(self) -> str:
        return f"range:{self.key}"


@dataclass(frozen=True)
class Rules:
    ranges: tuple[RangeRule, ...]
    required: tuple[str, ...]          # features that must be non-NaN

    @property
    def names(self) -> list[str]:
        """Every rule name in attribution order (for flow tables)."""
        return (["delineation"] + [f"missing:{f}" for f in (*self.required, *FIDUCIALS)]
                + [OVERLAP_RULE] + [r.name for r in self.ranges])


@dataclass
class QCResult:
    failed: list[tuple[str, str]] = field(default_factory=list)   # (rule, detail)

    @property
    def passed(self) -> bool:
        return not self.failed

    @property
    def rules(self) -> list[str]:
        return [r for r, _ in self.failed]


def build_rules(cfg: dict | None = None, registry: Registry = REGISTRY) -> Rules:
    """Range rules from config `qc:`, checked against the registry declarations.

    A plausibility range must contain the feature's declared normal range and its
    abnormal threshold; otherwise QC would exclude the very finding under study.
    Required features: every range-checked feature plus the primary endpoint."""
    import pulse_ecg.features  # noqa: F401  (registers every feature)

    cfg = cfg or config.load()
    ranges = []
    for key, (lo, hi) in cfg["qc"].items():
        if key not in QC_FEATURES:
            raise ValueError(f"qc.{key}: no registered feature mapped in qc.QC_FEATURES")
        spec = registry[QC_FEATURES[key]]
        lo, hi = float(lo), float(hi)
        if not lo < hi:
            raise ValueError(f"qc.{key}: range must satisfy lo < hi, got {[lo, hi]}")
        if spec.normal is not None and not (lo <= spec.normal[0] and spec.normal[1] <= hi):
            raise ValueError(f"qc.{key} {[lo, hi]} excludes part of {spec.name}'s "
                             f"normal range {list(spec.normal)}")
        if spec.abnormal is not None and not lo <= spec.abnormal.value <= hi:
            raise ValueError(f"qc.{key} {[lo, hi]} excludes {spec.name}'s abnormal "
                             f"threshold {spec.abnormal}")
        ranges.append(RangeRule(key, spec.name, lo, hi, spec.units))
    required = [r.feature for r in ranges]
    primary = cfg.get("analysis", {}).get("primary", {}).get("feature")
    if primary and primary not in required:
        registry[primary]                                  # KeyError if unregistered
        required.append(primary)
    return Rules(tuple(ranges), tuple(required))


def _nan(v) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v))


def evaluate(row: Mapping[str, Any], rules: Rules | None = None) -> QCResult:
    """Every rule `row` fails. `row` is an extract()-style row plus the median-beat
    fiducials p_offset and qrs_onset (see qc_row). Group status is never read."""
    rules = rules or build_rules()
    res = QCResult()
    if row.get("delineation_reason"):
        res.failed.append(("delineation", str(row["delineation_reason"])))

    for f in rules.required:
        if _nan(row.get(f)):
            res.failed.append((f"missing:{f}", row.get(f"{f}_reason") or "value is NaN"))
    for p in FIDUCIALS:
        if _nan(row.get(p)):
            res.failed.append((f"missing:{p}", f"{p} not delineated on the median beat"))

    # Atrial depolarisation is complete at the P offset, so ventricular depolarisation
    # (QRS onset) cannot begin before it: qrs_onset < p_offset is physiologically
    # impossible and means the delineation of the median beat is wrong. Measured on
    # 15.5 % of Design A records (paper/notes/qrs_onset_overlap.md); in those records
    # QRS duration is inflated, PR shortened and the P axis flips by ~180 degrees, so
    # the record fails. Equality (P offset = QRS onset, PR segment 0) is allowed.
    pof, qon = row.get("p_offset"), row.get("qrs_onset")
    if not (_nan(pof) or _nan(qon)) and qon < pof:
        res.failed.append((OVERLAP_RULE, f"qrs_onset {qon:g} < p_offset {pof:g} (samples)"))

    for r in rules.ranges:
        v = row.get(r.feature)
        if not _nan(v) and not r.lo <= v <= r.hi:
            res.failed.append((r.name, f"{r.feature} {v:g} {r.units} outside "
                                       f"[{r.lo:g}, {r.hi:g}]"))
    return res


def apply(rows: pd.DataFrame, rules: Rules | None = None) -> pd.DataFrame:
    """`rows` with qc_pass, qc_failed (rule names, ';'-joined), qc_detail and
    qc_first (the rule a flow diagram attributes the exclusion to). Same length
    and order as the input."""
    rules = rules or build_rules()
    order = {n: i for i, n in enumerate(rules.names)}
    out = rows.copy()
    results = [evaluate(r, rules) for r in rows.to_dict("records")]
    out["qc_pass"] = [r.passed for r in results]
    out["qc_failed"] = [";".join(r.rules) for r in results]
    out["qc_detail"] = [" | ".join(f"{n}: {d}" for n, d in r.failed) for r in results]
    out["qc_first"] = [min(r.rules, key=order.__getitem__) if r.failed else ""
                       for r in results]
    return out


def flow(qc: pd.DataFrame, rules: Rules | None = None, by: str | None = None) -> pd.DataFrame:
    """Exclusion counts per rule: `any` = records failing the rule (a record may fail
    several), `first` = records attributed to it in rule order (sums to the total
    excluded). Split by `by` (e.g. group) when given; `by` is only read here, after
    QC has been decided."""
    rules = rules or build_rules()
    groups = [("all", qc)] + (list(qc.groupby(by)) if by else [])
    rows = []
    for g, sub in groups:
        n = len(sub)
        failed = sub["qc_failed"].str.split(";")
        rows.append({"group": g, "rule": "records", "any": n, "first": n})
        for name in rules.names:
            rows.append({"group": g, "rule": name,
                         "any": int(failed.apply(lambda fs: name in fs).sum()),
                         "first": int((sub["qc_first"] == name).sum())})
        rows.append({"group": g, "rule": "excluded", "any": int((~sub["qc_pass"]).sum()),
                     "first": int((~sub["qc_pass"]).sum())})
        rows.append({"group": g, "rule": "final_n", "any": int(sub["qc_pass"].sum()),
                     "first": int(sub["qc_pass"].sum())})
    df = pd.DataFrame(rows)
    df["n"] = df["group"].map({g: len(s) for g, s in groups})
    for c in ("any", "first"):
        df[f"{c}_pct"] = 100 * df[c] / df["n"]
    return df


def qc_row(rec, cfg: dict | None = None) -> dict[str, Any]:
    """extract()-style feature row plus delineation status and the median-beat
    fiducials evaluate() needs."""
    import pulse_ecg.features  # noqa: F401
    from .extract import context

    cfg = cfg or config.load()
    ctx = context(rec, cfg)
    d = ctx.delineation
    row: dict[str, Any] = {"record_id": rec.record_id, "delineation_reason": d.reason}
    row.update({p: d.points.get(p, math.nan) for p in FIDUCIALS})
    row.update(REGISTRY.compute(ctx.rec, ctx))
    return row
