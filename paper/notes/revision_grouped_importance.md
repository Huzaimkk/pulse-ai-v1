# POST-HOC revision: grouped permutation importance (random forest)

**Post hoc.** `models.py` is itself post hoc, and this analysis was added at revision. No model was retrained differently, and no existing output was changed.

Script: `paper/notes/revision_grouped_importance.py`. Output: `out/revision_grouped_importance.csv`.

## Method

**Models.** The script calls `models.py`'s own `load()` and `nested("random_forest", …)`. That gives the same 775 records, feature sets A and B, RF grid, outer 5-fold StratifiedGroupKFold by patient (seed 42) and inner tuning. It asserts that the refit models reproduce the stored `out/models_metrics.csv` results:

| Set | Pooled AUROC | Fold-mean AUROC |
|---|---|---|
| A | 0.8302 | 0.8335 |
| B | 0.7456 | 0.7475 |

Both reproduced exactly. These are therefore the same fitted outer-fold models that `models_importance.csv` was computed from.

**Grouped permutation.** On each outer test fold, all columns of a group are permuted with **one shared row permutation**. This keeps the correlation within the group and stops a correlated partner from standing in for a permuted feature.

- Drop = fold AUROC − AUROC with the group permuted.
- Each fold uses 30 repeats (seed 42), as in `models.py`.
- The table reports the mean ± SD (ddof 1) across the 5 outer folds.

**Groups**

| Group | Members |
|---|---|
| QT-family | qt_ms, qtc_bazett_ms, qtc_fridericia_ms, tpe_qt_ratio |
| P-family (set A only) | p_duration_ms, ptfv1_uv_ms, p_axis_deg, iab_class, mvp_score |
| QRS-family | qrs_duration_ms, fqrs_present |
| Rate | mean_hr_bpm, rmssd_ms |
| Demographics | age, sex |

Everything else is permuted on its own: tpe_ms, n_beats and the 12 ST-J60 leads. Each feature appears in exactly one group (asserted). n_beats is closely tied to HR but was kept on its own as specified.

The last column, "Σ individual", is the sum of the members' single-feature importances from `out/models_importance.csv`, for comparison.

## Results: drop in AUROC (mean ± SD across 5 outer folds)

### Set A (all features), RF, AUROC 0.830

| Group | Δ AUROC | SD | Σ individual |
|---|---|---|---|
| P-family | **0.1632** | 0.0158 | 0.1667 |
| QT-family | **0.0290** | 0.0104 | 0.0161 |
| st_j60_v6_mv | 0.0135 | 0.0048 | 0.0140 |
| QRS-family | 0.0068 | 0.0061 | 0.0080 |
| st_j60_v5_mv | 0.0053 | 0.0061 | 0.0050 |
| Rate | 0.0032 | 0.0048 | 0.0037 |
| st_j60_v4_mv | 0.0021 | 0.0021 | 0.0030 |
| st_j60_v2_mv | 0.0013 | 0.0047 | 0.0010 |
| Demographics | 0.0012 | 0.0009 | 0.0014 |
| st_j60_v1_mv | 0.0009 | 0.0021 | 0.0006 |
| st_j60_avl_mv | 0.0008 | 0.0013 | 0.0004 |
| st_j60_i_mv | 0.0004 | 0.0025 | 0.0006 |
| st_j60_v3_mv | 0.0004 | 0.0029 | 0.0015 |
| n_beats | 0.0002 | 0.0008 | 0.0003 |
| st_j60_avr_mv | 0.0001 | 0.0013 | 0.0003 |
| st_j60_avf_mv | 0.0001 | 0.0017 | 0.0003 |
| tpe_ms | −0.0002 | 0.0016 | −0.0001 |
| st_j60_ii_mv | −0.0003 | 0.0009 | −0.0001 |
| st_j60_iii_mv | −0.0008 | 0.0021 | −0.0011 |

### Set B (no P-derived features), RF, AUROC 0.746

| Group | Δ AUROC | SD | Σ individual |
|---|---|---|---|
| QT-family | **0.0464** | 0.0189 | 0.0187 |
| st_j60_v6_mv | 0.0263 | 0.0057 | 0.0281 |
| QRS-family | 0.0088 | 0.0044 | 0.0100 |
| st_j60_v5_mv | 0.0065 | 0.0085 | 0.0072 |
| Rate | 0.0033 | 0.0033 | 0.0031 |
| st_j60_i_mv | 0.0029 | 0.0025 | 0.0032 |
| Demographics | 0.0027 | 0.0018 | 0.0033 |
| st_j60_v3_mv | 0.0024 | 0.0051 | 0.0031 |
| st_j60_avl_mv | 0.0022 | 0.0033 | 0.0016 |
| st_j60_v2_mv | 0.0022 | 0.0037 | 0.0026 |
| st_j60_v4_mv | 0.0015 | 0.0016 | 0.0018 |
| st_j60_avf_mv | 0.0012 | 0.0027 | 0.0009 |
| st_j60_ii_mv | 0.0008 | 0.0020 | 0.0013 |
| n_beats | 0.0004 | 0.0013 | 0.0002 |
| st_j60_avr_mv | 0.0002 | 0.0023 | −0.0001 |
| st_j60_v1_mv | −0.0001 | 0.0032 | −0.0000 |
| tpe_ms | −0.0008 | 0.0021 | −0.0013 |
| st_j60_iii_mv | −0.0009 | 0.0058 | −0.0011 |

## Summary (post hoc)

1. **Set A.** The P-family dominates (ΔAUROC 0.163 ± 0.016), and almost all of it is PTFV1. Shuffling PTFV1 alone gives 0.161 in `models_importance.csv`, and the group is barely larger. The ranking does not change when correlated features are grouped.
2. **QT-family.** Its single-feature importance was diluted by correlation. Grouped, it is about 1.8× (set A, 0.029 vs Σ 0.016) to 2.5× (set B, 0.046 vs Σ 0.019) the sum of its members' individual importances. It is clearly second in set A and first in set B. The other groups are roughly additive: grouped ≈ Σ individual.
3. **Small effects.** Rate, demographics and most single ST leads have drops within about 1–2 SD of zero across folds. Of the single ST leads, only V6 is clearly above zero in both sets (A 0.014 ± 0.005, B 0.026 ± 0.006). Age and sex carry almost no signal, as expected for matching variables.
