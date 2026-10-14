"""
POST-HOC REVISION SENSITIVITY ANALYSIS (not pre-specified in analysis_plan.md).

Design A, final analysed cohort (n=775 after QC). Counts LVH / bundle-branch /
fascicular / WPW / pacemaker statements per group, then reruns the
pre-specified OLS (marker ~ group + age + sex, group_comparison.ols_group) for
five markers under two variants:
  (a) excl  - drop records with LVH, CLBBB, CRBBB, WPW or pacemaker (both
              groups; no re-matching)
  (b) lvhcov - full cohort, LVH added as a covariate

Reads only; existing outputs are not touched.

    .venv/Scripts/python.exe paper/notes/revision_lvh_bbb.py

Inputs:  out/qc_designA.csv, out/cohort.csv, out/group_comparison_designA.csv,
         PTB-XL ptbxl_database.csv
Outputs: out/revision_lvh_bbb.csv, out/revision_lvh_bbb_prevalence.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import group_comparison as gc  # noqa: E402  (same code path as the original)
from pulse_ecg import io  # noqa: E402

CODES = ["LVH", "CLBBB", "CRBBB", "IRBBB", "ILBBB", "LAFB", "WPW"]
EXCLUDE = ["LVH", "CLBBB", "CRBBB", "WPW", "pacemaker"]
MARKERS = ["p_duration_ms", "ptfv1_uv_ms", "qrs_duration_ms",
           "qtc_fridericia_ms", "mean_hr_bpm"]


def ols_group_cov(df, feat, extra):
    """gc.ols_group with extra covariates appended; identical formula otherwise."""
    d = df[[feat, "case", "age", "sex", *extra]].dropna()
    X = np.column_stack([np.ones(len(d)), d.case, d.age, d.sex, *(d[e] for e in extra)])
    y = d[feat].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    se = np.sqrt(resid @ resid / dof * np.linalg.inv(X.T @ X)[1, 1])
    t = beta[1] / se
    tc = stats.t.ppf(0.975, dof)
    return beta[1], beta[1] - tc * se, beta[1] + tc * se, 2 * stats.t.sf(abs(t), dof)


def load():
    qc = pd.read_csv(ROOT / "out" / "qc_designA.csv")
    kept = qc[qc["qc_pass"]].drop(columns="group")
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort.design == "A"][["ecg_id", "group", "age", "sex"]]
    df = kept.merge(a, on="ecg_id", validate="one_to_one")
    df["case"] = (df.group == "case").astype(int)
    assert len(df) == 775 and df.case.sum() == 259, (len(df), df.case.sum())

    db = io.ptbxl_database()
    sub = db.loc[df.ecg_id]
    for c in CODES:                       # any stated likelihood, as in 01_cohort.py
        df[c] = sub.scp_codes.apply(lambda s, c=c: c in s).to_numpy().astype(int)
    df["pacemaker_col"] = sub.pacemaker.notna().to_numpy().astype(int)
    df["PACE_code"] = sub.scp_codes.apply(lambda s: "PACE" in s).to_numpy().astype(int)
    df["pacemaker"] = (df.pacemaker_col | df.PACE_code).astype(int)
    df["any_excluded"] = df[EXCLUDE].max(axis=1)
    return df


def prevalence(df):
    rows = []
    for c in [*CODES, "pacemaker_col", "PACE_code", "pacemaker", "any_excluded"]:
        r = {"item": c}
        for lab, g in (("case", 1), ("control", 0)):
            x = df.loc[df.case == g, c]
            r |= {f"{lab}_k": int(x.sum()), f"{lab}_n": len(x), f"{lab}_pct": 100 * x.mean()}
        rows.append(r)
    return pd.DataFrame(rows)


def fit(df, feat, extra=()):
    ca = df.loc[df.case == 1, feat].dropna()
    co = df.loc[df.case == 0, feat].dropna()
    diff, lo, hi, p = ols_group_cov(df, feat, list(extra)) if extra else gc.ols_group(df, feat)
    g = gc.hedges_g(ca, co)[0]
    return {"n_case": len(ca), "n_control": len(co), "adj_diff": diff,
            "ci_lo": lo, "ci_hi": hi, "p": p, "hedges_g": g}


def main():
    df = load()
    prev = prevalence(df)

    orig = pd.read_csv(ROOT / "out" / "group_comparison_designA.csv").set_index("feature")
    excl = df[df.any_excluded == 0]
    rows = []
    for f in MARKERS:
        o = fit(df, f)
        # guard: recomputed original must equal the stored result
        so = orig.loc[f]
        assert np.isclose(o["adj_diff"], so.adj_diff) and np.isclose(o["p"], so.p), f
        assert np.isclose(*(ols_group_cov(df, f, [])[0], o["adj_diff"]))
        for variant, r in (("original", o),
                           ("a_excl_LVH_CLBBB_CRBBB_WPW_pacemaker", fit(excl, f)),
                           ("b_full_LVH_covariate", fit(df, f, ["LVH"]))):
            rows.append({"analysis": "POST-HOC revision sensitivity",
                         "marker": f, "variant": variant, **r})
    res = pd.DataFrame(rows)

    out = ROOT / "out"
    res.to_csv(out / "revision_lvh_bbb.csv", index=False)
    prev.to_csv(out / "revision_lvh_bbb_prevalence.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print("POST-HOC revision sensitivity analysis (not pre-specified)\n")
        print(prev.round(1).to_string(index=False))
        print(f"\nexcluded in (a): {len(df) - len(excl)} records "
              f"({(df.any_excluded & df.case).sum()} cases, "
              f"{(df.any_excluded & (1 - df.case)).sum()} controls)\n")
        print(res.drop(columns="analysis").to_string(index=False,
              float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
