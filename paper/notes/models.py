"""
Stage 8 (POST-HOC, not in analysis_plan.md): binary classification, LAE case vs
control, Design A, QC-passed records (n = 775). The plan's "Not done" section says
no classifier is fitted; this script is a later addition and is labelled post-hoc.

    .venv/Scripts/python.exe paper/notes/models.py            # full run
    .venv/Scripts/python.exe paper/notes/models.py --n-perm 5 # quick check

Feature sets: A = all registered single-value features + age + sex;
              B = A minus the P-derived features (p_duration_ms, ptfv1_uv_ms,
                  iab_class, mvp_score, p_axis_deg).
The 300 s HRV features are excluded (NaN on 10 s records).

Models: L2 logistic regression (standardised), decision tree (max_depth 2-5),
random forest, HistGradientBoosting. Nested CV: outer 5-fold
StratifiedGroupKFold by patient_id (seed 42), inner 3-fold (same class, seed 42).
Median imputation, scaling, hyperparameter search and the decision threshold
(Youden J on inner out-of-fold probabilities) are all fitted inside the outer
training fold only.

Inputs:  out/qc_designA.csv, out/cohort.csv
Outputs: out/models_metrics.csv, out/models_confusion.csv,
         out/models_importance.csv, out/models_permutation_test.csv,
         out/models_oof_predictions.csv, out/decision_tree_set{A,B}.txt,
         figures/fig/fig_decision_tree.pdf (set A), fig_decision_tree_setB.pdf,
         paper/notes/models.md
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, confusion_matrix, roc_auc_score,
                             roc_curve)
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, export_text, plot_tree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "figures"))
import pulse_ecg.features  # noqa: E402,F401  (registers every feature)
from pulse_ecg.registry import REGISTRY  # noqa: E402

SEED = 42
N_BOOT = 1000
OUT = ROOT / "out"
FIG = ROOT / "figures" / "fig"

LONG_WINDOW = ["sdnn_ms", "pnn50_pct", "lf_ms2", "hf_ms2", "lf_hf_ratio", "sampen",
               "dfa_alpha1"]                    # require >= 300 s; NaN on 10 s records
NOT_SCALAR = ["rr_intervals_ms"]                # a list per record, not a feature value
CIRCULAR = ["p_duration_ms", "ptfv1_uv_ms", "iab_class", "mvp_score", "p_axis_deg"]
DEMOG = ["age", "sex"]

MODELS = {
    "logistic_l2": (
        Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler()),
                  ("clf", LogisticRegression(penalty="l2", max_iter=5000))]),
        {"clf__C": [0.01, 0.1, 1.0, 10.0]}),
    "decision_tree": (
        Pipeline([("imp", SimpleImputer(strategy="median")),
                  ("clf", DecisionTreeClassifier(random_state=SEED))]),
        {"clf__max_depth": [2, 3, 4, 5]}),
    "random_forest": (
        Pipeline([("imp", SimpleImputer(strategy="median")),
                  ("clf", RandomForestClassifier(n_estimators=300, random_state=SEED))]),
        {"clf__max_features": ["sqrt", 0.5], "clf__min_samples_leaf": [1, 5, 20]}),
    "hist_gradient_boosting": (
        Pipeline([("imp", SimpleImputer(strategy="median")),
                  ("clf", HistGradientBoostingClassifier(random_state=SEED))]),
        {"clf__learning_rate": [0.03, 0.1], "clf__max_depth": [2, 3]}),
}


# ------------------------------------------------------------------- data ---
def load():
    qc = pd.read_csv(OUT / "qc_designA.csv")
    kept = qc[qc["qc_pass"]].drop(columns="group")
    cohort = pd.read_csv(OUT / "cohort.csv")
    a = cohort[cohort.design == "A"][["ecg_id", "group", "patient_id", "age", "sex"]]
    df = kept.merge(a, on="ecg_id", validate="one_to_one")
    df["case"] = (df.group == "case").astype(int)
    assert len(df) == 775 and df.case.sum() == 259, (len(df), df.case.sum())
    assert df[LONG_WINDOW].isna().all().all()

    reg = [s.name for s in REGISTRY]
    feats_a = [f for f in reg if f not in LONG_WINDOW + NOT_SCALAR] + DEMOG
    assert set(CIRCULAR) <= set(feats_a)
    feats_b = [f for f in feats_a if f not in CIRCULAR]
    return df, {"A": feats_a, "B": feats_b}


# ------------------------------------------------------------- nested CV ---
def outer_cv(y, g, seed=SEED):
    return list(StratifiedGroupKFold(5, shuffle=True, random_state=seed)
                .split(np.zeros(len(y)), y, g))


def inner_cv():
    return StratifiedGroupKFold(3, shuffle=True, random_state=SEED)


def fit_fold(name, X, y, g, tr, n_jobs=1):
    """Tune on the outer training fold, pick the threshold there, refit there."""
    est, grid = MODELS[name]
    gs = GridSearchCV(clone(est), grid, scoring="roc_auc", cv=inner_cv(), n_jobs=n_jobs,
                      refit=True)
    gs.fit(X[tr], y[tr], groups=g[tr])
    best = clone(est).set_params(**gs.best_params_)
    oof = cross_val_predict(best, X[tr], y[tr], groups=g[tr], cv=inner_cv(),
                            method="predict_proba", n_jobs=n_jobs)[:, 1]
    fpr, tpr, thr = roc_curve(y[tr], oof)
    j = np.argmax(tpr - fpr)
    thresh = min(thr[j], 1.0)                 # thr[0] is inf; never selected in practice
    return gs.best_estimator_, gs.best_params_, thresh


def nested(name, X, y, g, n_jobs=1, seed=SEED):
    prob = np.full(len(y), np.nan)
    pred = np.full(len(y), -1)
    fold = np.full(len(y), -1)
    info = []
    for k, (tr, te) in enumerate(outer_cv(y, g, seed)):
        model, params, thresh = fit_fold(name, X, y, g, tr, n_jobs)
        prob[te] = model.predict_proba(X[te])[:, 1]
        pred[te] = (prob[te] >= thresh).astype(int)
        fold[te] = k
        info.append({"fold": k, "params": params, "threshold": thresh, "model": model,
                     "te": te, "auroc": roc_auc_score(y[te], prob[te])})
    assert not np.isnan(prob).any() and (pred >= 0).all()
    return prob, pred, fold, info


# ---------------------------------------------------------------- metrics ---
def threshold_metrics(y, pred):
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    prec = tp / (tp + fp) if tp + fp else np.nan
    f1 = 2 * tp / (2 * tp + fp + fn)
    return {"accuracy": (tp + tn) / len(y), "sensitivity": sens, "specificity": spec,
            "precision": prec, "f1": f1, "tn": tn, "fp": fp, "fn": fn, "tp": tp}


def bootstrap_ci(y, p, rng):
    """Percentile CI. Every patient has exactly one record, so resampling records
    is resampling patients."""
    n = len(y)
    au, ap = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        if y[i].min() == y[i].max():
            continue
        au.append(roc_auc_score(y[i], p[i]))
        ap.append(average_precision_score(y[i], p[i]))
    assert len(au) >= 0.99 * N_BOOT
    return np.percentile(au, [2.5, 97.5]), np.percentile(ap, [2.5, 97.5])


def check_row(r):
    """Guards against impossible numbers before anything is written."""
    for k in ("accuracy", "sensitivity", "specificity", "f1", "auroc", "auprc"):
        assert 0 <= r[k] <= 1, (r["model"], k, r[k])
    prec = 0.0 if np.isnan(r["precision"]) else r["precision"]
    assert r["f1"] <= max(prec, r["sensitivity"]) + 1e-12, \
        (r["model"], r["set"], r["f1"], prec, r["sensitivity"])
    assert r["auroc_lo"] <= r["auroc"] <= r["auroc_hi"]
    assert r["auprc_lo"] <= r["auprc"] <= r["auprc_hi"]
    assert r["tn"] + r["fp"] + r["fn"] + r["tp"] == 775


def baseline_row(y, fset):
    prev = y.mean()
    pred = np.zeros_like(y)                           # always predict the majority (control)
    r = {"set": fset, "model": "majority_class", "auroc": 0.5, "auroc_lo": 0.5,
         "auroc_hi": 0.5, "auprc": prev, "auprc_lo": prev, "auprc_hi": prev,
         "fold_auroc_mean": 0.5, "fold_auroc_sd": 0.0, "threshold_median": np.nan,
         "best_params": "", **threshold_metrics(y, pred)}
    return r


# ----------------------------------------------------- permutation tests ---
def perm_once(name, X, y, g, seed):
    yp = np.random.default_rng(seed).permutation(y)
    prob, *_ = nested(name, X, yp, g)
    return roc_auc_score(yp, prob)


def perm_importance(info, X, y, feats):
    rows = []
    for f in info:
        te = f["te"]
        pi = permutation_importance(f["model"], X[te], y[te], scoring="roc_auc",
                                    n_repeats=30, random_state=SEED, n_jobs=-1)
        rows.append(pi.importances_mean)
    m = np.array(rows)
    return (pd.DataFrame({"feature": feats, "importance_mean": m.mean(0),
                          "importance_sd_across_folds": m.std(0, ddof=1)})
            .sort_values("importance_mean", ascending=False, ignore_index=True))


# ------------------------------------------------------------------- tree ---
def full_tree(X, y, g, feats, fset):
    est, grid = MODELS["decision_tree"]
    cv = StratifiedGroupKFold(5, shuffle=True, random_state=SEED)
    gs = GridSearchCV(clone(est), grid, scoring="roc_auc", cv=cv, n_jobs=-1).fit(X, y, groups=g)
    tree = gs.best_estimator_.named_steps["clf"]
    txt = (f"Decision tree, feature set {fset}, refit on all 775 records "
           f"(max_depth {gs.best_params_['clf__max_depth']} chosen by 5-fold grouped CV "
           f"AUROC {gs.best_score_:.3f}). Missing values median-imputed.\n"
           f"class 0 = control, class 1 = LAE case\n\n"
           + export_text(tree, feature_names=feats, show_weights=True))
    (OUT / f"decision_tree_set{fset}.txt").write_text(txt, encoding="utf-8")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import figstyle as S
    S.apply()
    leaves, depth = tree.get_n_leaves(), tree.get_depth()
    fig, ax = plt.subplots(figsize=(max(S.COL_2, 1.6 * leaves), 1.25 * (depth + 1)))
    plot_tree(tree, feature_names=feats, class_names=["control", "LAE"], filled=False,
              impurity=False, proportion=False, rounded=True, fontsize=8, ax=ax)   # >= S.MIN_PT; must be int
    FIG.mkdir(parents=True, exist_ok=True)
    stem = "fig_decision_tree" if fset == "A" else "fig_decision_tree_setB"
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    return gs.best_params_["clf__max_depth"], gs.best_score_, depth, leaves


# ---------------------------------------------------------------- report ---
def fmt(v, d=3):
    return "n/a" if pd.isna(v) else f"{v:.{d}f}"


def write_md(met, conf, imp, perm, trees, feats, n_perm):
    L = ["# Stage 8: classification models, LAE case vs control (Design A)", "",
         "**POST-HOC. Not pre-specified.** `analysis_plan.md` (\"Not done\") states that "
         "no predictive model, classifier or decision tree is fitted. These models were "
         "added afterwards at the author's request and must be reported as post-hoc "
         "and exploratory.", "",
         "Generated by `paper/notes/models.py`. Do not edit numbers by hand.", "",
         f"n = 775 QC-passed records (259 cases, 516 controls; prevalence "
         f"{259/775:.3f}), one record per patient. Nested CV: outer 5-fold "
         "StratifiedGroupKFold by `patient_id` (seed 42), inner 3-fold for "
         "hyperparameters and threshold. Median imputation and scaling are fitted inside "
         "each training fold. Threshold = Youden J on inner out-of-fold probabilities. "
         "AUROC/AUPRC are computed on the pooled outer-fold predictions, with 95% "
         "percentile bootstrap CIs (1000 resamples). The fold-mean AUROC ± SD is also shown.", "",
         f"- Set A ({len(feats['A'])} features): all registered single-value features "
         "+ age + sex.",
         f"- Set B ({len(feats['B'])} features): set A without p_duration_ms, "
         "ptfv1_uv_ms, iab_class, mvp_score, p_axis_deg (the P-derived features; the "
         "LAE label was read from these same P waves).",
         "- Age and sex are matching variables (1:2, sex exact, age ±5 y), so they "
         "carry almost no signal by design.",
         "- fqrs_present is missing on 290 of 775 records and is median-imputed.", "",
         "## Results (outer folds)", "",
         "| Set | Model | AUROC (95% CI) | AUPRC (95% CI) | Fold AUROC mean ± SD | "
         "Accuracy | Sensitivity | Specificity | Precision | F1 |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in met.iterrows():
        L.append(f"| {r['set']} | {r['model']} | {fmt(r.auroc)} ({fmt(r.auroc_lo)}–"
                 f"{fmt(r.auroc_hi)}) | {fmt(r.auprc)} ({fmt(r.auprc_lo)}–{fmt(r.auprc_hi)}) "
                 f"| {fmt(r.fold_auroc_mean)} ± {fmt(r.fold_auroc_sd)} | {fmt(r.accuracy)} | "
                 f"{fmt(r.sensitivity)} | {fmt(r.specificity)} | {fmt(r.precision)} | "
                 f"{fmt(r.f1)} |")
    L += ["", "Majority-class baseline: always predicts control; AUROC 0.5, AUPRC = "
          "prevalence. Its precision is undefined (no positive predictions).", ""]

    L += ["## Pooled confusion matrix, best model per set", ""]
    for _, p in perm.iterrows():
        c = conf[(conf.set == p.set) & (conf.model == p.model)].iloc[0]
        L += [f"Set {p.set}, {p.model} (threshold from the inner loop, pooled over 5 "
              "outer folds):", "",
              "| | Predicted control | Predicted LAE |", "|---|---|---|",
              f"| Control | {c.tn} | {c.fp} |", f"| LAE | {c.fn} | {c.tp} |", ""]

    L += ["## Permutation test of AUROC", "",
          "| Set | Model | Observed AUROC | Permutations | Mean null AUROC | p |",
          "|---|---|---|---|---|---|"]
    for _, p in perm.iterrows():
        L.append(f"| {p.set} | {p.model} | {fmt(p.observed_auroc)} | {p.n_perm} | "
                 f"{fmt(p.null_mean)} | {fmt(p.p_value, 4)} |")
    L += ["", f"Labels permuted {n_perm} times; the whole nested CV (tuning and "
          "threshold) is rerun on each. p = (1 + #null ≥ observed) / (1 + n). The model "
          "tested was chosen as the best of four, so this p is not adjusted for that "
          "selection.", ""]

    L += ["## Permutation feature importance, best model per set (top 10)", "",
          "Drop in AUROC when a feature is shuffled, measured on each outer test fold "
          "with that fold's fitted model (30 repeats), mean over the 5 folds.", ""]
    for _, p in perm.iterrows():
        t = imp[(imp.set == p.set)].head(10)
        L += [f"Set {p.set}, {p.model}:", "", "| Rank | Feature | ΔAUROC mean | SD across folds |",
              "|---|---|---|---|"]
        for i, r in enumerate(t.itertuples(), 1):
            L.append(f"| {i} | {r.feature} | {r.importance_mean:.4f} | "
                     f"{r.importance_sd_across_folds:.4f} |")
        L.append("")

    L += ["## Decision tree refit on all 775", ""]
    for fset, (md, sc, depth, leaves) in trees.items():
        stem = "fig_decision_tree" if fset == "A" else "fig_decision_tree_setB"
        L.append(f"- Set {fset}: max_depth {md} (5-fold grouped CV AUROC {sc:.3f}); "
                 f"depth {depth}, {leaves} leaves. `figures/fig/{stem}.pdf`, "
                 f"`out/decision_tree_set{fset}.txt`.")
    L.append("")

    ba = met[(met.set == "A") & (met.model != "majority_class")].sort_values("auroc").iloc[-1]
    bb = met[(met.set == "B") & (met.model != "majority_class")].sort_values("auroc").iloc[-1]
    L += ["## Set A vs set B", "",
          f"Best AUROC is {ba.auroc:.3f} ({ba.auroc_lo:.3f}–{ba.auroc_hi:.3f}, {ba.model}) "
          f"with all features and {bb.auroc:.3f} ({bb.auroc_lo:.3f}–{bb.auroc_hi:.3f}, "
          f"{bb.model}) without the P-derived features "
          f"(difference {ba.auroc - bb.auroc:+.3f}; the CIs are from separate bootstraps and "
          "are not a paired test).", ""]
    (ROOT / "paper" / "notes" / "models.md").write_text("\n".join(L), encoding="utf-8")


# ------------------------------------------------------------------ main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=200)
    args = ap.parse_args()

    df, feats = load()
    y = df.case.to_numpy()
    g = df.patient_id.to_numpy()
    assert abs(y.mean() - 0.334) < 0.001

    met_rows, conf_rows, oof_rows, best_info = [], [], [], {}
    for fset, fl in feats.items():
        X = df[fl].to_numpy(float)
        met_rows.append(baseline_row(y, fset))
        for name in MODELS:
            prob, pred, fold, info = nested(name, X, y, g, n_jobs=-1)
            (alo, ahi), (plo, phi) = bootstrap_ci(y, prob, np.random.default_rng(SEED))
            fa = [f["auroc"] for f in info]
            r = {"set": fset, "model": name, "auroc": roc_auc_score(y, prob),
                 "auroc_lo": alo, "auroc_hi": ahi,
                 "auprc": average_precision_score(y, prob), "auprc_lo": plo, "auprc_hi": phi,
                 "fold_auroc_mean": np.mean(fa), "fold_auroc_sd": np.std(fa, ddof=1),
                 "threshold_median": np.median([f["threshold"] for f in info]),
                 "best_params": "; ".join(str(f["params"]) for f in info),
                 **threshold_metrics(y, pred)}
            met_rows.append(r)
            best_info[(fset, name)] = info
            oof_rows.append(pd.DataFrame({"ecg_id": df.ecg_id, "set": fset, "model": name,
                                          "fold": fold, "case": y, "prob": prob, "pred": pred}))
            print(f"{fset} {name:24s} AUROC {r['auroc']:.3f} ({alo:.3f}-{ahi:.3f}) "
                  f"AUPRC {r['auprc']:.3f} F1 {r['f1']:.3f}", flush=True)

    met = pd.DataFrame(met_rows)
    for _, r in met.iterrows():
        check_row(r)
    conf = met[["set", "model", "tn", "fp", "fn", "tp"]]

    # Model-stage outputs are written before the (slow) permutation tests, so an
    # interrupted run still leaves them on disk.
    best_model, imp_rows = {}, []
    for fset, fl in feats.items():
        X = df[fl].to_numpy(float)
        sub = met[(met.set == fset) & (met.model != "majority_class")]
        best_model[fset] = best = sub.sort_values("auroc").iloc[-1]
        imp = perm_importance(best_info[(fset, best.model)], X, y, fl)
        imp_rows.append(imp.assign(set=fset, model=best.model))
    imp = pd.concat(imp_rows, ignore_index=True)

    met.to_csv(OUT / "models_metrics.csv", index=False)
    conf.to_csv(OUT / "models_confusion.csv", index=False)
    imp.to_csv(OUT / "models_importance.csv", index=False)
    pd.concat(oof_rows, ignore_index=True).to_csv(OUT / "models_oof_predictions.csv",
                                                  index=False)
    print("wrote models_metrics/confusion/importance/oof_predictions.csv", flush=True)

    perm_rows = []
    perm_path = OUT / "models_permutation_test.csv"
    perm_path.unlink(missing_ok=True)                 # no stale rows from an earlier run
    for fset, fl in feats.items():
        X = df[fl].to_numpy(float)
        best = best_model[fset]
        null = []
        for v in Parallel(n_jobs=-1, return_as="generator")(       # ordered, same seeds
                delayed(perm_once)(best.model, X, y, g, SEED + 1 + i)
                for i in range(args.n_perm)):
            null.append(v)
            if len(null) % 20 == 0:
                print(f"{fset} permutation {len(null)}/{args.n_perm}", flush=True)
        null = np.array(null)
        perm_rows.append({"set": fset, "model": best.model, "observed_auroc": best.auroc,
                          "n_perm": args.n_perm, "null_mean": null.mean(),
                          "null_sd": null.std(ddof=1), "null_max": null.max(),
                          "p_value": (1 + (null >= best.auroc).sum()) / (1 + args.n_perm)})
        pd.DataFrame(perm_rows).to_csv(perm_path, index=False)
        print(f"{fset} permutation test {best.model}: p = {perm_rows[-1]['p_value']:.4f}",
              flush=True)
    perm = pd.DataFrame(perm_rows)

    trees = {fset: full_tree(df[fl].to_numpy(float), y, g, fl, fset)
             for fset, fl in feats.items()}

    write_md(met, conf, imp, perm, trees, feats, args.n_perm)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(met.drop(columns="best_params").round(3).to_string())
        print(perm.round(4).to_string())


if __name__ == "__main__":
    main()
