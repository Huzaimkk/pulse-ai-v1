# POST-HOC revision sensitivity analysis: LVH, bundle-branch block, WPW, pacemaker

**This analysis is post hoc. It was added at revision and is not pre-specified in `analysis_plan.md`.** It does not replace the pre-specified Design A result (`out/group_comparison_designA.csv`), and neither that file nor the plan was changed.

Script: `paper/notes/revision_lvh_bbb.py`. Outputs: `out/revision_lvh_bbb.csv` (models) and `out/revision_lvh_bbb_prevalence.csv` (counts).
Cohort: the final analysed Design A cohort, n = 775 after QC (259 cases, 516 controls).

## 1. Conduction and structural statements in the analysed cohort

A statement counts as present when its SCP code appears in `scp_codes` at any likelihood. That is the same rule `01_cohort.py` uses (`min_likelihood: 0`).

| Statement | Cases n (%) of 259 | Controls n (%) of 516 |
|---|---|---|
| LVH | 67 (25.9) | 44 (8.5) |
| CLBBB | 26 (10.0) | 17 (3.3) |
| CRBBB | 19 (7.3) | 17 (3.3) |
| IRBBB | 13 (5.0) | 22 (4.3) |
| ILBBB | 3 (1.2) | 1 (0.2) |
| LAFB | 30 (11.6) | 58 (11.2) |
| WPW | 0 (0.0) | 1 (0.2) |
| Pacemaker present* | 0 (0.0) | 1 (0.2) |
| Any of LVH/CLBBB/CRBBB/WPW/pacemaker | 111 (42.9) | 78 (15.1) |

\*"Pacemaker present" means a non-empty `pacemaker` column or a `PACE` SCP code. No record in the cohort has the `PACE` code. The one control flagged (ecg_id 1625) has the column value `PACE????, nan`, which is uncertain. That record also has CRBBB, so variant (a) excludes it either way.

**Were these exclusion criteria?** No. `scripts/01_cohort.py` does not exclude any of these statements. Its only selection rules are a sinus code (`db["sinus"] = db.scp_codes.apply(has_code, args=(c["sinus_codes"], 0))`), the label or never-label definition, one record per patient, file present, age ≥ 40 (`usable = lambda d: d[d.file_present & (d.age >= c["min_age"])]`) and sex/age matching. The terms LVH, BBB, LAFB, WPW and pacemaker do not appear in `scripts/`, `config.yaml`, `pulse_ecg/` or `paper/notes/qc_exclusions.py`.

## 2. Pre-specified OLS model (marker ~ group + age + sex) under two post-hoc variants

- **Original**: the pre-specified result. The script recomputes it and asserts that it equals the stored `group_comparison_designA.csv`.
- **(a) Exclusion**: records with LVH, CLBBB, CRBBB, WPW or pacemaker are removed from both groups (111 cases, 78 controls). Pairs are not re-matched.
- **(b) LVH covariate**: the full cohort, with LVH (0/1) added to the model.

The difference is case minus control, adjusted for age and sex (and LVH in b). Hedges' g is the unadjusted between-group effect size, computed as in the original code, so in (b) it equals the original g by construction.

| Marker | Variant | n case | n ctrl | Adj. diff | 95% CI | p | Hedges' g |
|---|---|---|---|---|---|---|---|
| P duration (ms), primary | original | 259 | 516 | 2.04 | −0.88 to 4.95 | 0.170 | 0.101 |
| | (a) exclusion | 148 | 438 | 3.91 | 0.37 to 7.46 | 0.031 | 0.188 |
| | (b) +LVH | 259 | 516 | 1.58 | −1.42 to 4.58 | 0.301 | 0.101 |
| PTFV1 (µV·ms) | original | 258 | 516 | −3627 | −4237 to −3017 | 4.3e−29 | −0.882 |
| | (a) exclusion | 148 | 438 | −3430 | −4153 to −2706 | 2.5e−19 | −0.877 |
| | (b) +LVH | 258 | 516 | −3351 | −3974 to −2728 | 1.9e−24 | −0.882 |
| QRS (ms) | original | 259 | 516 | 14.16 | 10.37 to 17.94 | 5.5e−13 | 0.544 |
| | (a) exclusion | 148 | 438 | 9.52 | 5.64 to 13.41 | 1.9e−06 | 0.436 |
| | (b) +LVH | 259 | 516 | 16.60 | 12.77 to 20.43 | 8.8e−17 | 0.544 |
| QTcF (ms) | original | 259 | 516 | 32.75 | 24.95 to 40.56 | 7.7e−16 | 0.619 |
| | (a) exclusion | 148 | 438 | 32.38 | 22.80 to 41.96 | 7.3e−11 | 0.615 |
| | (b) +LVH | 259 | 516 | 33.62 | 25.59 to 41.66 | 8.8e−16 | 0.619 |
| HR (bpm) | original | 259 | 516 | 4.94 | 3.15 to 6.74 | 8.2e−08 | 0.411 |
| | (a) exclusion | 148 | 438 | 3.64 | 1.44 to 5.85 | 0.0012 | 0.307 |
| | (b) +LVH | 259 | 516 | 5.00 | 3.16 to 6.85 | 1.3e−07 | 0.411 |

The p values are unadjusted for multiplicity and come from a post-hoc analysis.

## Summary (post hoc)

1. No marker changed direction in either variant. PTFV1, QRS, QTcF and HR stayed significant, and their effect sizes changed little; the largest change was QRS in (a), from 14.2 to 9.5 ms.
2. The primary P-duration result lost no significance, because it was not significant to begin with (p = 0.170). It stayed non-significant with LVH as a covariate (p = 0.301). In the exclusion subset it became nominally significant (+3.9 ms, p = 0.031). That subset is post hoc, unmatched and uncorrected for multiplicity, so it does not change the pre-specified null primary result.
3. LVH and left/right BBB were 2–3× more common in cases than in controls (any of them: 42.9% vs 15.1%). Even so, the secondary group differences hold after these records are excluded or LVH is adjusted for.
