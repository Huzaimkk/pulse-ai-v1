"""
POST-HOC REVISION ANALYSIS (not pre-specified; models.py is itself post hoc).

Grouped permutation importance for the random forest, feature sets A and B.
The models are refit with models.py's own nested() (same RF, grid, outer
StratifiedGroupKFold folds, seed 42). The script asserts that each fold AUROC
and the pooled AUROC equal the stored out/models_metrics.csv values, so these
are the same fitted models. Existing outputs are not rewritten.

Correlated features are permuted TOGETHER: one row permutation, applied to
every column of the group, so within-group correlation is kept and the group
cannot "cover" for itself. Drop = fold AUROC - AUROC with the group permuted;
30 repeats per fold (as models.py), mean per fold, then mean +/- SD (ddof 1)
across the 5 outer folds.

    .venv/Scripts/python.exe paper/notes/revision_grouped_importance.py

Output: out/revision_grouped_importance.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import models as M  # noqa: E402  (same data, folds and estimators)

N_REPEATS = 30
GROUPS = {
    "QT-family": ["qt_ms", "qtc_bazett_ms", "qtc_fridericia_ms", "tpe_qt_ratio"],
    "P-family": ["p_duration_ms", "ptfv1_uv_ms", "p_axis_deg", "iab_class", "mvp_score"],
    "QRS-family": ["qrs_duration_ms", "fqrs_present"],
    "rate": ["mean_hr_bpm", "rmssd_ms"],
    "demographics": ["age", "sex"],
}


def groups_for(feats):
    out, used = {}, set()
    for name, members in GROUPS.items():
        m = [f for f in members if f in feats]
        if m:
            out[name] = m
            used |= set(m)
    for f in feats:
        if f not in used:
            out[f] = [f]
    assert sorted(sum(out.values(), [])) == sorted(feats)   # every feature exactly once
    return out


def grouped_drop(model, X, y, cols, rng):
    base = roc_auc_score(y, model.predict_proba(X)[:, 1])
    drops = []
    for _ in range(N_REPEATS):
        Xp = X.copy()
        Xp[:, cols] = X[rng.permutation(len(X))][:, cols]
        drops.append(base - roc_auc_score(y, model.predict_proba(Xp)[:, 1]))
    return float(np.mean(drops))


def main():
    df, feats = M.load()
    y, g = df.case.to_numpy(), df.patient_id.to_numpy()
    met = pd.read_csv(ROOT / "out" / "models_metrics.csv")
    ind = pd.read_csv(ROOT / "out" / "models_importance.csv")
    rows = []
    for fset, fl in feats.items():
        X = df[fl].to_numpy(float)
        prob, _, _, info = M.nested("random_forest", X, y, g, n_jobs=-1)
        stored = met[(met.set == fset) & (met.model == "random_forest")].iloc[0]
        assert np.isclose(roc_auc_score(y, prob), stored.auroc), fset
        assert np.isclose(np.mean([f["auroc"] for f in info]), stored.fold_auroc_mean), fset
        grp = groups_for(fl)
        idx = {f: i for i, f in enumerate(fl)}
        per_fold = {name: [] for name in grp}
        for f in info:
            te = f["te"]
            rng = np.random.default_rng(M.SEED)
            for name, members in grp.items():
                per_fold[name].append(grouped_drop(f["model"], X[te], y[te],
                                                   [idx[m] for m in members], rng))
        iset = ind[ind.set == fset].set_index("feature")
        for name, members in grp.items():
            v = np.array(per_fold[name])
            rows.append({"analysis": "POST-HOC revision", "set": fset, "model": "random_forest",
                         "group": name, "members": "; ".join(members),
                         "n_members": len(members),
                         "drop_auroc_mean": v.mean(), "drop_auroc_sd_across_folds": v.std(ddof=1),
                         **{f"fold{k}": x for k, x in enumerate(v)},
                         "sum_individual_importance (models_importance.csv)":
                             iset.loc[members, "importance_mean"].sum()})
        print(f"set {fset}: RF reproduced (pooled AUROC {stored.auroc:.4f}, "
              f"fold mean {stored.fold_auroc_mean:.4f})", flush=True)

    res = (pd.DataFrame(rows)
           .sort_values(["set", "drop_auroc_mean"], ascending=[True, False], ignore_index=True))
    res.to_csv(ROOT / "out" / "revision_grouped_importance.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20,
                           "display.max_colwidth", 70):
        print(res[["set", "group", "members", "drop_auroc_mean", "drop_auroc_sd_across_folds",
                   "sum_individual_importance (models_importance.csv)"]]
              .to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
