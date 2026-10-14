"""
make_numbers.py — every number quoted in the paper text, as LaTeX macros.

    python make_numbers.py --out-dir ../out > numbers.tex

Rule 8: no number enters the paper unless a re-runnable script produced it. The
.tex file uses only the macros written here (\\nPrimDiff, \\nAucA, ...), so changing
an upstream result and re-running this script updates the paper.

Values not in a CSV are taken from the notes that the analysis scripts wrote; each
is marked SOURCE below.
"""

import argparse
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", required=True, type=Path)
o = ap.parse_args().out_dir
M = {}


def f(x, d=1):
    return f"{x:.{d}f}"


def sgn(x, d=1):
    s = f"{x:+.{d}f}"
    return s.replace("-", "$-$") if s.startswith("-") else s


def p_fmt(p):
    if p < 0.001:
        return "$<$0.001"
    return f"{p:.3f}" if p < 0.01 else f"{p:.2f}"


# ------------------------------------------------------------------ cohort ---
g1 = pd.read_csv(o / "gate1.csv").set_index("design")
M["nEligCase"] = f"{int(g1.loc['A', 'cases_eligible']):,}"
M["nEligCtrl"] = f"{int(g1.loc['A', 'controls_eligible']):,}"
M["nMatchCase"] = int(g1.loc["A", "cases_matched"])
M["nMatchCtrl"] = int(g1.loc["A", "controls_matched"])
M["nMatchAll"] = M["nMatchCase"] + M["nMatchCtrl"]
M["ageMatched"] = f(g1.loc["A", "mean_age_matched"])
M["nBCaseMatched"] = int(g1.loc["B", "cases_matched"])
M["nBCtrlMatched"] = int(g1.loc["B", "controls_matched"])

fl = pd.read_csv(o / "qc_flow_designA.csv")
q = lambda g, r, c="any": fl[(fl.group == g) & (fl.rule == r)][c].iloc[0]
M["nExcl"] = int(q("all", "excluded"))
M["pctExcl"] = f(q("all", "excluded", "any_pct"))
M["nOverlap"] = int(q("all", "qrs_onset_before_p_offset"))
M["pctOverlap"] = f(q("all", "qrs_onset_before_p_offset", "any_pct"))
M["pctOverlapCase"] = f(q("case", "qrs_onset_before_p_offset", "any_pct"))
M["pctOverlapCtrl"] = f(q("control", "qrs_onset_before_p_offset", "any_pct"))
M["nQrsRange"] = int(q("all", "range:qrs_ms", "first"))
M["nPRange"] = int(q("all", "range:p_ms", "first"))
M["nQtRange"] = int(q("all", "range:qt_ms", "first"))
M["nCase"] = int(q("case", "final_n"))
M["nCtrl"] = int(q("control", "final_n"))
M["nAll"] = M["nCase"] + M["nCtrl"]
M["prev"] = f(M["nCase"] / M["nAll"], 3)

t1 = pd.read_csv(o / "table1_postqc_designA.csv")
post = t1[t1.stage == "post-QC"].set_index("group")
M["ageCase"] = f(post.loc["case", "age_mean"]); M["ageCaseSd"] = f(post.loc["case", "age_sd"])
M["ageCtrl"] = f(post.loc["control", "age_mean"]); M["ageCtrlSd"] = f(post.loc["control", "age_sd"])
M["maleCase"] = f(post.loc["case", "male_pct"]); M["maleCtrl"] = f(post.loc["control", "male_pct"])
M["smdAge"] = f(post.loc["SMD", "age_mean"], 3)
M["smdMale"] = sgn(post.loc["SMD", "male_pct"], 3)

pv = pd.read_csv(o / "posthoc_report_provenance_designA.csv")
pq = pv[pv.set == "QC-passed"].set_index("group")
M["nHumanCase"] = int(pq.loc["case", "validated_by_human"])
M["nHumanCtrl"] = int(pq.loc["control", "validated_by_human"])

# --------------------------------------------------------------- gate 2 ---
s = pd.read_csv(o / "gate2_summary.csv")
lt = s[(s.dataset == "ludb") & (s.half == "test") & (s.reference == "knn") & (s.rhythm == "sinus")]
lt = lt.set_index("boundary")
for b, k in [("p_onset", "Pon"), ("p_offset", "Poff"), ("qrs_onset", "Qon"),
             ("qrs_offset", "Qoff"), ("t_offset", "Toff")]:
    M[f"sd{k}"] = f(lt.loc[b, "sd_within_record_ms"])
    M[f"mean{k}"] = sgn(lt.loc[b, "mean_ms"])
    M[f"mart{k}"] = f(lt.loc[b, "martinez2004_sd_ms"])
    M[f"cse{k}"] = f(lt.loc[b, "cse_2sd_tol_ms"])
M["nLudbBeats"] = int(lt.loc["p_onset", "n"])
M["pSens"] = f(100 * lt.loc["p_onset", "p_sensitivity"])

pdur = pd.read_csv(o / "gate2_p_duration.csv")
r = pdur[(pdur.dataset == "ludb") & (pdur.half == "test") & (pdur.reference == "knn")].iloc[0]
M["nLudbRec"] = int(r.n_records)
M["gTwoBias"] = sgn(r.record_mean_err_ms); M["gTwoBiasSd"] = f(r.record_sd_err_ms)
M["gTwoMae"] = f(r.record_mae_ms)
M["gTwoRefPrev"] = f(100 * r.ref_prev_ge_thr, 0); M["gTwoGotPrev"] = f(100 * r.got_prev_ge_thr, 0)
M["gTwoKappa"] = f(r.kappa_ge_thr, 2)
r = pdur[pdur.dataset == "qtdb"].iloc[0]
M["qtdbRec"] = int(r.n_records)
M["qtdbBias"] = sgn(r.record_mean_err_ms); M["qtdbBiasSd"] = f(r.record_sd_err_ms)
M["qtdbMae"] = f(r.record_mae_ms)

b = pd.read_csv(o / "gate2_beats.csv")
t = b[(b.dataset == "ludb") & (b.half == "test") & (b.reference == "knn") & (b.sinus == 1)]
t = t.dropna(subset=["p_onset_err_ms", "p_offset_err_ms", "p_dur_ref_ms"])
e = t.p_offset_err_ms - t.p_onset_err_ms
M["baBias"] = sgn(e.mean()); M["baLo"] = sgn(e.mean() - 1.96 * e.std(ddof=1), 0)
M["baHi"] = sgn(e.mean() + 1.96 * e.std(ddof=1), 0)

# SOURCE: paper/notes/p_duration_compression.md (written by its measurement script)
M["compSlope"] = "0.782"; M["compIcpt"] = "27.9"

# --------------------------------------------------------------- gate 3 ---
g3 = pd.read_csv(o / "gate3_offsets_designA.csv")
g3 = g3.set_index(["feature", "stratum"])
M["offPLo"] = f(g3.loc[("p_duration_ms", "pooled"), "offset_lo"], 0)
M["offPHi"] = f(g3.loc[("p_duration_ms", "pooled"), "offset_hi"], 0)
M["offQLo"] = f(g3.loc[("qrs_duration_ms", "pooled"), "offset_lo"], 0)
M["offQHi"] = f(g3.loc[("qrs_duration_ms", "pooled"), "offset_hi"], 0)
M["medP"] = f(g3.loc[("p_duration_ms", "pooled"), "median"], 0)
M["medQ"] = f(g3.loc[("qrs_duration_ms", "pooled"), "median"], 0)

# ---------------------------------------------------------- group comparison ---
g = pd.read_csv(o / "group_comparison_designA.csv").set_index("feature")
P = g.loc["p_duration_ms"]
M["pdCaseMed"] = f(P.case_median, 0); M["pdCtrlMed"] = f(P.control_median, 0)
M["primDiff"] = sgn(P.adj_diff); M["primLo"] = sgn(P.adj_diff_lo); M["primHi"] = sgn(P.adj_diff_hi)
M["primP"] = p_fmt(P.p); M["primG"] = sgn(P.hedges_g, 2)
M["primGlo"] = sgn(P.g_lo, 2); M["primGhi"] = sgn(P.g_hi, 2); M["primCliff"] = f(P.cliffs_delta, 2)
for feat, k, d in [("ptfv1_uv_ms", "Ptf", 0), ("qrs_duration_ms", "Qrs", 1),
                   ("qtc_fridericia_ms", "QtcF", 1), ("qtc_bazett_ms", "QtcB", 1),
                   ("mean_hr_bpm", "Hr", 1), ("tpe_ms", "Tpe", 1), ("p_axis_deg", "Pax", 1),
                   ("st_j60_v6_mv", "StVsix", 3), ("st_j60_avr_mv", "StAvr", 3),
                   ("rmssd_ms", "Rmssd", 1)]:
    r = g.loc[feat]
    M[f"diff{k}"] = sgn(r.adj_diff, d); M[f"lo{k}"] = sgn(r.adj_diff_lo, d)
    M[f"hi{k}"] = sgn(r.adj_diff_hi, d); M[f"g{k}"] = sgn(r.hedges_g, 2)
    M[f"q{k}"] = p_fmt(r.q_bh)
sec = g[g.role != "primary"].dropna(subset=["q_bh"])
M["nSec"] = len(sec); M["nSecSig"] = int((sec.q_bh < 0.05).sum())

th = pd.read_csv(o / "group_threshold_calls_designA.csv")
th = th.set_index(["feature", "call"])
M["iabCase"] = f(th.loc[("iab_class", "any IAB (>=1)"), "case_pct"])
M["iabCtrl"] = f(th.loc[("iab_class", "any IAB (>=1)"), "control_pct"])
M["ptfAbCase"] = f(th.loc[("ptfv1_uv_ms", "abnormal (<= -4000)"), "case_pct"])
M["ptfAbCtrl"] = f(th.loc[("ptfv1_uv_ms", "abnormal (<= -4000)"), "control_pct"])
M["mvpHiCase"] = f(th.loc[("mvp_score", "high (>= 5)"), "case_pct"])
M["mvpHiCtrl"] = f(th.loc[("mvp_score", "high (>= 5)"), "control_pct"])

# ------------------------------------------------------------- sensitivity ---
s1 = pd.read_csv(o / "sens1_matched_sets_designA.csv")
c = [x for x in s1.columns]
row = s1[s1.feature == "p_duration_ms"].iloc[0]
M["sOneN"] = int(row.n_sets)
M["sOneDiff"] = sgn(row.mean_d); M["sOneLo"] = sgn(row.d_lo)
M["sOneHi"] = sgn(row.d_hi); M["sOneP"] = p_fmt(row.p_paired_t)

gb = pd.read_csv(o / "group_comparison_designB.csv").set_index("feature")
B = gb.loc["p_duration_ms"]
M["bCase"] = int(B.n_case); M["bCtrl"] = int(B.n_control)
M["bDiff"] = sgn(B.adj_diff); M["bLo"] = sgn(B.adj_diff_lo); M["bHi"] = sgn(B.adj_diff_hi)
M["bP"] = p_fmt(B.p); M["bG"] = sgn(B.hedges_g, 2)
B2 = gb.loc["ptfv1_uv_ms"]
M["bPtfDiff"] = sgn(B2.adj_diff, 0); M["bPtfQ"] = p_fmt(B2.q_bh)

# SOURCE: paper/notes/sensitivity.md section 3 (post-hoc, human-validated subset)
M["hvDiff"] = "+1.8"; M["hvLo"] = "$-$1.6"; M["hvHi"] = "+5.3"; M["hvP"] = "0.29"

# ------------------------------------------------------------------ models ---
m = pd.read_csv(o / "models_metrics.csv").set_index(["set", "model"])
for st in "AB":
    for mdl, k in [("logistic_l2", "Lr"), ("decision_tree", "Dt"),
                   ("random_forest", "Rf"), ("hist_gradient_boosting", "Gb")]:
        r = m.loc[(st, mdl)]
        p_ = f"{k}{st}"
        M[f"auc{p_}"] = f(r.auroc, 3); M[f"aucLo{p_}"] = f(r.auroc_lo, 3)
        M[f"aucHi{p_}"] = f(r.auroc_hi, 3); M[f"prc{p_}"] = f(r.auprc, 3)
        M[f"sens{p_}"] = f(r.sensitivity, 3); M[f"spec{p_}"] = f(r.specificity, 3)
        M[f"acc{p_}"] = f(r.accuracy, 3); M[f"fone{p_}"] = f(r.f1, 3)
        M[f"prec{p_}"] = f(r.precision, 3)
        assert r.f1 <= max(r.precision, r.sensitivity) + 1e-9, (st, mdl)
    M[f"accMaj{st}"] = f(m.loc[(st, "majority_class")].accuracy, 3)
r = m.loc[("A", "random_forest")]
M["tn"], M["fp"], M["fn"], M["tp"] = (int(r.tn), int(r.fp), int(r.fn), int(r.tp))
pt = pd.read_csv(o / "models_permutation_test.csv").set_index("set")
M["nPerm"] = int(pt.loc["A", "n_perm"])
M["permP"] = f(pt.p_value.max(), 3)
M["nullMax"] = f(pt.null_max.max(), 2)
M["aucGapAB"] = f(m.loc[("A", "random_forest")].auroc - m.loc[("B", "random_forest")].auroc, 3)


# ------------------------------------------------- POST-HOC revision (Oct 2026) ---
# Every value below is post hoc and is labelled so in the paper.
lv = pd.read_csv(o / "revision_lvh_bbb.csv")
lv = lv.set_index(["marker", "variant"])
EX, CV = "a_excl_LVH_CLBBB_CRBBB_WPW_pacemaker", "b_full_LVH_covariate"
r = lv.loc[("p_duration_ms", EX)]
M["exN"] = int(r.n_case); M["exNCtrl"] = int(r.n_control)
for feat, k, d in [("p_duration_ms", "P", 1), ("ptfv1_uv_ms", "Ptf", 0),
                   ("qrs_duration_ms", "Qrs", 1), ("qtc_fridericia_ms", "QtcF", 1),
                   ("mean_hr_bpm", "Hr", 1)]:
    for v, tag in [(EX, "ex"), (CV, "cov")]:
        r = lv.loc[(feat, v)]
        M[f"{tag}{k}Diff"] = sgn(r.adj_diff, d); M[f"{tag}{k}Lo"] = sgn(r.ci_lo, d)
        M[f"{tag}{k}Hi"] = sgn(r.ci_hi, d); M[f"{tag}{k}P"] = p_fmt(r.p)

pr = pd.read_csv(o / "revision_lvh_bbb_prevalence.csv").set_index("item")
for it, k in [("LVH", "Lvh"), ("CLBBB", "Lbbb"), ("CRBBB", "Rbbb")]:
    M[f"prev{k}Case"] = f(pr.loc[it, "case_pct"]); M[f"prev{k}Ctrl"] = f(pr.loc[it, "control_pct"])
M["exDropCase"] = int(pr.loc["any_excluded", "case_k"])
M["exDropCtrl"] = int(pr.loc["any_excluded", "control_k"])

qi = pd.read_csv(o / "revision_qc_included_excluded.csv")
kc = qi[(qi.group == "key check") & (qi.variable == "Exclusion rate")].iloc[0]
M["qcRateP"] = p_fmt(kc.p)
fc = pd.read_csv(o / "revision_qc_flow_chain.csv")
lost = lambda g, st: int(fc[(fc.group == g) & (fc.step == st)].lost.iloc[0])
M["ageLostCase"] = lost("case", "lost: age < 40")
M["qcLostCase"] = lost("case", "lost: QC exclusion")
M["qcLostCtrl"] = lost("control", "lost: QC exclusion")
M["qcRateCase"] = f(100 * M["qcLostCase"] / M["nMatchCase"])
M["qcRateCtrl"] = f(100 * M["qcLostCtrl"] / M["nMatchCtrl"])

pl = pd.read_csv(o / "revision_ptfv1_ludb.csv")
H = 'Gate2 held-out half, sinus'
gp = lambda sub, met, lvl="per-record": pl[(pl.subset == sub) & (pl.metric == met) & (pl.level == lvl)].iloc[0]
r = gp(H, "P offset error (ours - LUDB), ms")
M["vOffBias"] = sgn(r.bias); M["vOffSd"] = f(r.sd); M["vOffN"] = int(r.n)
r = gp(H, "P onset error (ours - LUDB), ms")
M["vOnBias"] = sgn(r.bias); M["vOnSd"] = f(r.sd)
r = gp(H, "PTFV1 ours - LUDB-boundaries, µV·ms")
M["vPtfBias"] = sgn(r.bias, 0); M["vPtfLo"] = sgn(r.loa_lo, 0); M["vPtfHi"] = sgn(r.loa_hi, 0)
M["vIcc"] = f(gp(H, "ICC(2,1)").value, 2); M["vRho"] = f(gp(H, "Spearman r").value, 2)
r = gp(H, "Cohen kappa, abnormal <= -4000")
M["vKappa"] = f(r.value, 2)
M["vOursOnly"] = int(r.note.split("ours only ")[1].split(",")[0])
M["vLudbOnly"] = int(r.note.split("LUDB only ")[1])
r = gp(H, "terminal negative duration, ms ours - LUDB"); M["vTermDur"] = sgn(r.bias)
r = gp(H, "terminal negative depth, µV ours - LUDB"); M["vTermDepth"] = sgn(r.bias)
M["vAbnOurs"] = f(gp(H, "PTFV1 <= -4000, ours").pct, 0)
M["vAbnLudb"] = f(gp(H, "PTFV1 <= -4000, LUDB boundaries").pct, 0)
C = "PTB-XL analysed controls"
M["ctlTermDepth"] = f(gp(C, "terminal depth median, µV").value, 0)
M["ludbTermDepth"] = f(gp("all sinus", "terminal depth median, µV (LUDB, ours)").value, 0)

gi = pd.read_csv(o / "revision_grouped_importance.csv").set_index(["set", "group"])
for st in "AB":
    r = gi.loc[(st, "QT-family")]
    M[f"giQt{st}"] = f(r.drop_auroc_mean, 3); M[f"giQtSum{st}"] = f(r.iloc[-1], 3)
r = gi.loc[("A", "P-family")]
M["giPA"] = f(r.drop_auroc_mean, 3); M["giPSdA"] = f(r.drop_auroc_sd_across_folds, 3)

# SOURCE: paper/notes/revision_methods_facts.md (post-hoc measurements and code facts)
M["bBefore"] = 51; M["bAfter"] = 67
M["notchBpDb"] = "17.9"; M["notchMedDb"] = "0.5"; M["notchOverTen"] = 6
M["qtdbBeats"] = "3,178"; M["qtdbScored"] = 104; M["qtdbListed"] = 105

# ------------------------------------------------------------------ write ---
import re
print("% numbers.tex — generated by make_numbers.py. Do not edit by hand.")
for k, v in M.items():
    assert re.fullmatch(r"[A-Za-z]+", k), k
    print(f"\\newcommand{{\\n{k[0].upper() + k[1:]}}}{{{v}}}")
