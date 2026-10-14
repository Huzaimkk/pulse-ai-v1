"""
POST-HOC REVISION ANALYSIS (not pre-specified in analysis_plan.md).

Design A: records that entered QC (981 matched records), INCLUDED (qc_pass) vs
EXCLUDED, overall and by group, plus the eligible -> usable -> matched ->
analysed reconciliation. QC rules and existing outputs are not touched; QC
results are read from out/qc_designA.csv as written by qc_exclusions.py.

    .venv/Scripts/python.exe paper/notes/revision_qc_included_excluded.py

Outputs: out/revision_qc_included_excluded.csv (characteristics),
         out/revision_qc_reasons.csv, out/revision_qc_flow_chain.csv
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pulse_ecg import config, io  # noqa: E402

RULE_ORDER = ["delineation", "missing:mean_hr_bpm", "missing:qrs_duration_ms",
              "missing:qt_ms", "missing:p_duration_ms", "missing:p_offset",
              "missing:qrs_onset", "qrs_onset_before_p_offset", "range:hr",
              "range:qrs_ms", "range:qt_ms", "range:p_ms"]
CONT = [("age", "Age, y"), ("mean_hr_bpm", "HR, bpm"), ("n_beats", "Beats, n"),
        ("p_duration_ms", "P duration, ms"), ("qrs_duration_ms", "QRS, ms"),
        ("qtc_fridericia_ms", "QTcF, ms")]


def load_cohort_module():
    spec = importlib.util.spec_from_file_location("cohort01", ROOT / "scripts" / "01_cohort.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def chi2(table):
    """Pearson chi-square without continuity correction; Fisher exact alongside."""
    table = np.asarray(table)
    chi, p, *_ = stats.chi2_contingency(table, correction=False)
    pf = stats.fisher_exact(table)[1] if table.shape == (2, 2) else np.nan
    return chi, p, pf


def characteristics(df):
    rows = []
    for grp in ("all", "case", "control"):
        d = df if grp == "all" else df[df.group == grp]
        inc, exc = d[d.qc_pass], d[~d.qc_pass]
        rows.append({"group": grp, "variable": "n", "included": f"{len(inc)}",
                     "excluded": f"{len(exc)}", "test": "", "p": np.nan})
        for col, lab in CONT:
            a, b = inc[col].dropna(), exc[col].dropna()
            p = stats.ttest_ind(a, b, equal_var=False).pvalue
            rows.append({"group": grp, "variable": lab,
                         "included": f"{a.mean():.1f} ± {a.std():.1f} (n={len(a)})",
                         "excluded": f"{b.mean():.1f} ± {b.std():.1f} (n={len(b)})",
                         "test": "Welch t", "p": p})
        ki, ke = int(inc.male.sum()), int(exc.male.sum())
        _, p, _ = chi2([[ki, len(inc) - ki], [ke, len(exc) - ke]])
        rows.append({"group": grp, "variable": "Male, n (%)",
                     "included": f"{ki} ({100 * ki / len(inc):.1f})",
                     "excluded": f"{ke} ({100 * ke / len(exc):.1f})",
                     "test": "chi-square", "p": p})
    return pd.DataFrame(rows)


def reasons(df):
    rows = []
    for rule in RULE_ORDER:
        r = {"rule": rule}
        for grp in ("all", "case", "control"):
            d = df if grp == "all" else df[df.group == grp]
            r[f"{grp}_first"] = int((d.qc_first == rule).sum())
            r[f"{grp}_any"] = int(d.qc_failed.fillna("").str.contains(rule, regex=False).sum())
        rows.append(r)
    out = pd.DataFrame(rows)
    tot = {"rule": "TOTAL excluded"}
    for grp in ("all", "case", "control"):
        d = df if grp == "all" else df[df.group == grp]
        n_exc = int((~d.qc_pass).sum())
        assert out[f"{grp}_first"].sum() == n_exc, grp      # flow reasons sum exactly
        tot[f"{grp}_first"] = n_exc
        tot[f"{grp}_any"] = np.nan
    return pd.concat([out, pd.DataFrame([tot])], ignore_index=True)


def chain(df):
    """332 eligible -> usable -> matched -> analysed, with the reason for each loss."""
    m = load_cohort_module()
    c = config.load()["cohort"]
    db = m.eligible(io.ptbxl_database(), c, c["min_likelihood"])
    cases, controls = m.design_a(db)
    rows = []
    for grp, d in (("case", cases), ("control", controls)):
        no_file = ~d.file_present
        young = d.file_present & (d.age < c["min_age"])
        usable = d.file_present & (d.age >= c["min_age"])
        matched = int((df.group == grp).sum())
        analysed = int(((df.group == grp) & df.qc_pass).sum())
        rows += [
            {"group": grp, "step": "eligible (sinus, first record per patient)", "n": len(d), "lost": np.nan},
            {"group": grp, "step": "lost: records500 file absent", "n": np.nan, "lost": int(no_file.sum())},
            {"group": grp, "step": f"lost: age < {c['min_age']}", "n": np.nan, "lost": int(young.sum())},
            {"group": grp, "step": "usable", "n": int(usable.sum()), "lost": np.nan},
            {"group": grp, "step": "lost: not selected by matching", "n": np.nan,
             "lost": int(usable.sum()) - matched},
            {"group": grp, "step": "matched (entered QC)", "n": matched, "lost": np.nan},
            {"group": grp, "step": "lost: QC exclusion", "n": np.nan, "lost": matched - analysed},
            {"group": grp, "step": "analysed", "n": analysed, "lost": np.nan},
        ]
        assert len(d) - no_file.sum() - young.sum() == usable.sum()
    return pd.DataFrame(rows)


def main():
    qc = pd.read_csv(ROOT / "out" / "qc_designA.csv")
    cohort = pd.read_csv(ROOT / "out" / "cohort.csv")
    a = cohort[cohort.design == "A"][["ecg_id", "age", "sex"]]
    df = qc.merge(a, on="ecg_id", validate="one_to_one")
    df["male"] = (df.sex == 0).astype(int)                 # PTB-XL: 0 = male
    assert len(df) == 981 and df.qc_pass.sum() == 775

    char = characteristics(df)
    rs = reasons(df)
    ch = chain(df)

    # key check: exclusion rate case vs control
    ec, nc = int((~df.qc_pass & (df.group == "case")).sum()), int((df.group == "case").sum())
    ek, nk = int((~df.qc_pass & (df.group == "control")).sum()), int((df.group == "control").sum())
    chi, p, pf = chi2([[ec, nc - ec], [ek, nk - ek]])
    # composition of exclusions (first reason) case vs control
    comp = rs.iloc[:-1].set_index("rule")[["case_first", "control_first"]]
    comp = comp[comp.sum(axis=1) > 0]
    chi_c, p_c, _ = chi2(comp.to_numpy())
    key = pd.DataFrame([
        {"group": "key check", "variable": "Exclusion rate",
         "included": f"case {ec}/{nc} ({100 * ec / nc:.1f}%)",
         "excluded": f"control {ek}/{nk} ({100 * ek / nk:.1f}%)",
         "test": f"chi-square={chi:.3f} (Fisher p={pf:.3g})", "p": p},
        {"group": "key check", "variable": "First-reason mix (4 reasons x group)",
         "included": "case " + ", ".join(f"{k}={v}" for k, v in comp.case_first.items()),
         "excluded": "control " + ", ".join(f"{k}={v}" for k, v in comp.control_first.items()),
         "test": f"chi-square={chi_c:.3f}, df={(len(comp) - 1)}", "p": p_c}])
    out = pd.concat([char, key], ignore_index=True)
    out.insert(0, "analysis", "POST-HOC revision")

    out.to_csv(ROOT / "out" / "revision_qc_included_excluded.csv", index=False)
    rs.to_csv(ROOT / "out" / "revision_qc_reasons.csv", index=False)
    ch.to_csv(ROOT / "out" / "revision_qc_flow_chain.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20,
                           "display.max_colwidth", 120):
        print("POST-HOC revision: QC included vs excluded, Design A\n")
        print(out.drop(columns="analysis").to_string(index=False,
              float_format=lambda v: f"{v:.3g}"))
        print("\n", rs.to_string(index=False))
        print("\n", ch.to_string(index=False))


if __name__ == "__main__":
    main()
