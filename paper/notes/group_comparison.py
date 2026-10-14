"""
Case-control comparison, Design A, exactly as pre-specified in
paper/notes/analysis_plan.md. No model, classifier or sensitivity analysis.

    .venv/Scripts/python.exe paper/notes/group_comparison.py

Inputs:  out/qc_designA.csv (qc_exclusions.py), out/cohort.csv
Outputs: out/group_comparison_designA.csv, out/group_threshold_calls_designA.csv
"""

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pulse_ecg.features  # noqa: E402,F401  (registers every feature)
from pulse_ecg.registry import REGISTRY  # noqa: E402

PRIMARY = "p_duration_ms"
ST = [f"st_j60_{l}_mv" for l in
      ("i", "ii", "iii", "avr", "avl", "avf", "v1", "v2", "v3", "v4", "v5", "v6")]
SECONDARY = ["n_beats", "mean_hr_bpm", "rmssd_ms", "qrs_duration_ms", "qt_ms",
             "qtc_bazett_ms", "qtc_fridericia_ms", "tpe_ms", "tpe_qt_ratio", *ST,
             "ptfv1_uv_ms", "p_axis_deg"]
# Not endpoints (plan: data-independent reasons)
LONG_WINDOW = ["sdnn_ms", "pnn50_pct", "lf_ms2", "hf_ms2", "lf_hf_ratio", "sampen",
               "dfa_alpha1"]
DESCRIPTIVE = ["iab_class", "mvp_score", "fqrs_present", "rr_intervals_ms"]

CIRCULAR = {"p_duration_ms": "circular", "ptfv1_uv_ms": "circular",
            "iab_class": "circular", "mvp_score": "circular",
            "p_axis_deg": "P-derived, not LAE criterion"}

PTFV1_ABN, MVP_HIGH = -4000, 5          # config.yaml atrial thresholds


def ols_group(df, feat):
    """feature ~ group + age + sex; classical SEs, t on n-4 df."""
    d = df[[feat, "case", "age", "sex"]].dropna()
    X = np.column_stack([np.ones(len(d)), d.case, d.age, d.sex])
    y = d[feat].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    se = np.sqrt(resid @ resid / dof * np.linalg.inv(X.T @ X)[1, 1])
    t = beta[1] / se
    tc = stats.t.ppf(0.975, dof)
    return beta[1], beta[1] - tc * se, beta[1] + tc * se, 2 * stats.t.sf(abs(t), dof)


def hedges_g(a, b):
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    g = (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (n1 + n2) - 9))
    se = np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))
    return g, g - 1.959964 * se, g + 1.959964 * se


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    b = np.sort(b)
    gt = np.searchsorted(b, a, side="left").sum()         # pairs with b < a
    lt = len(b) * len(a) - np.searchsorted(b, a, side="right").sum()   # b > a
    return (gt - lt) / (len(a) * len(b))


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1)
    return out


def row(df, feat, role):
    ca = df.loc[df.case == 1, feat].dropna()
    co = df.loc[df.case == 0, feat].dropna()
    diff, lo, hi, p = ols_group(df, feat)
    g, glo, ghi = hedges_g(ca, co)
    return {"feature": feat, "role": role,
            "circularity": CIRCULAR.get(feat, "not circular"),
            "n_case": len(ca), "n_control": len(co),
            "case_median": ca.median(), "case_q1": ca.quantile(.25), "case_q3": ca.quantile(.75),
            "control_median": co.median(), "control_q1": co.quantile(.25),
            "control_q3": co.quantile(.75),
            "adj_diff": diff, "adj_diff_lo": lo, "adj_diff_hi": hi, "p": p,
            "hedges_g": g, "g_lo": glo, "g_hi": ghi, "cliffs_delta": cliffs_delta(ca, co)}


def pct(x, cond):
    x = x.dropna()
    return int(cond(x).sum()), len(x), 100 * cond(x).mean()


def threshold_calls(df):
    calls = [
        ("iab_class", "none (0)", lambda v: v == 0),
        ("iab_class", "partial (1)", lambda v: v == 1),
        ("iab_class", "advanced (2)", lambda v: v == 2),
        ("iab_class", "any IAB (>=1)", lambda v: v >= 1),
        ("ptfv1_uv_ms", f"abnormal (<= {PTFV1_ABN})", lambda v: v <= PTFV1_ABN),
        ("mvp_score", "low (0-2)", lambda v: v <= 2),
        ("mvp_score", "intermediate (3-4)", lambda v: v.between(3, 4)),
        ("mvp_score", f"high (>= {MVP_HIGH})", lambda v: v >= MVP_HIGH),
        ("fqrs_present", "present (1)", lambda v: v == 1),
    ]
    rows = []
    for feat, call, cond in calls:
        r = {"feature": feat, "call": call, "circularity": CIRCULAR.get(feat, "not circular")}
        for lab, c in (("case", 1), ("control", 0)):
            k, n, pc = pct(df.loc[df.case == c, feat], cond)
            r |= {f"{lab}_k": k, f"{lab}_n": n, f"{lab}_pct": pc}
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    plan = ROOT / "paper" / "notes" / "analysis_plan.md"
    print("analysis_plan.md sha256:", hashlib.sha256(plan.read_bytes()).hexdigest())

    # Every registered feature is assigned exactly once.
    assigned = [PRIMARY, *SECONDARY, *LONG_WINDOW, *DESCRIPTIVE]
    assert sorted(assigned) == sorted(s.name for s in REGISTRY), \
        set(s.name for s in REGISTRY) ^ set(assigned)

    qc = pd.read_csv(ROOT / "out" / "qc_designA.csv")
    kept = qc[qc["qc_pass"]].drop(columns="group")
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort.design == "A"][["ecg_id", "group", "age", "sex"]]
    df = kept.merge(a, on="ecg_id", validate="one_to_one")
    df["case"] = (df.group == "case").astype(int)
    assert len(df) == 775 and df.case.sum() == 259, (len(df), df.case.sum())
    assert df[LONG_WINDOW].isna().all().all()

    res = pd.DataFrame([row(df, PRIMARY, "primary")] +
                       [row(df, f, "secondary") for f in SECONDARY])
    sec = res.role == "secondary"
    res["q_bh"] = np.nan
    res.loc[sec, "q_bh"] = bh(res.loc[sec, "p"])

    # mvp_score: descriptive distribution only (no model, no effect size)
    mv = {"feature": "mvp_score", "role": "descriptive", "circularity": "circular"}
    for lab, c in (("case", 1), ("control", 0)):
        v = df.loc[df.case == c, "mvp_score"].dropna()
        mv |= {f"n_{lab}": len(v), f"{lab}_median": v.median(),
               f"{lab}_q1": v.quantile(.25), f"{lab}_q3": v.quantile(.75)}
    res = pd.concat([res, pd.DataFrame([mv])], ignore_index=True)

    res.to_csv(ROOT / "out" / "group_comparison_designA.csv", index=False)
    calls = threshold_calls(df)
    calls.to_csv(ROOT / "out" / "group_threshold_calls_designA.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 99):
        print(res.round(4).to_string())
        print(calls.round(2).to_string())


if __name__ == "__main__":
    main()
