# Pre-specified analysis plan: case-control comparison, Design A

Written 2026-10-10. This plan was written and saved before any group statistic
(per-group summary, difference, test or effect size) was computed for any feature.
Before writing it, the only things inspected were: the post-QC counts and Table I
(`out/table1_postqc_designA.csv`), the matched-set structure after QC (counts of
cases and controls per `match_id`), pooled missingness counts per feature, the
pooled value counts of `iab_class`, `mvp_score` and `fqrs_present`, and the
feature registry and `config.yaml`. None of these were split by group.

The plan text below stays as written once Step B (running the analysis) starts. Anything that turns
out to be impossible or that has to change goes in **Deviations** at the end. The
plan text above that section is never edited.

## Data

- Design A, QC-passed records only: `out/qc_designA.csv` with `qc_pass == True`,
  n = 775 (259 cases, 516 controls). Age and sex come from `out/cohort.csv`
  (sex coded 0 = male, 1 = female, as in PTB-XL).
- Cases: sinus-rhythm PTB-XL records labelled LAO/LAE. Controls: sinus-rhythm
  records from patients who never carry LAO/LAE. Each case was matched 1:2 on sex
  (exact) and age (±5 y).
- Missing values: each feature is analysed on the records where it is non-missing
  (complete case per feature). n per group is reported for every feature. Values
  are not imputed.

## Primary endpoint

- `p_duration_ms`: global multilead P-wave duration, treated as continuous.
- Primary model: ordinary least-squares linear regression,
  `p_duration_ms ~ group + age + sex`, where group = 1 for case and 0 for control,
  age is linear in years, and sex is an indicator. Classical OLS standard errors
  are used, with a t distribution on n − 4 df.
- Reported estimate: the adjusted mean difference (cases minus controls) in ms,
  with a two-sided 95% CI and p value.
- Why unmatched: QC broke matched sets (159 complete 1:2 sets remain; the rest are
  1:1, case-only or control-only). The two groups are therefore analysed unmatched,
  with age and sex as covariates. After QC they remain balanced on both: Table I
  post-QC gives SMD 0.006 for age and −0.028 for male sex.
- The primary endpoint is tested on its own at α = 0.05, two-sided, with no
  multiplicity correction.

## Effect sizes (every continuous endpoint)

- Hedges' g: (mean_case − mean_control) / pooled SD, multiplied by the small-sample
  correction J = 1 − 3 / (4(n1 + n2) − 9). The 95% CI uses the large-sample normal
  approximation, SE(g) = sqrt((n1 + n2)/(n1·n2) + g² / (2(n1 + n2))). This value is
  unadjusted.
- Cliff's delta: P(case > control) − P(case < control) over all case-control pairs.
  It is unadjusted and robust to skew. Only the point estimate is reported.
- Median (IQR) per group, which is descriptive.

## Secondary endpoints

All registered continuous features other than `p_duration_ms`, each analysed with
the same model (`feature ~ group + age + sex`) and the same effect sizes:

| feature | family |
|---|---|
| n_beats | rate |
| mean_hr_bpm | rate |
| rmssd_ms | rate |
| qrs_duration_ms | ventricular |
| qt_ms | ventricular |
| qtc_bazett_ms | ventricular |
| qtc_fridericia_ms | ventricular |
| tpe_ms | ventricular |
| tpe_qt_ratio | ventricular |
| st_j60_{i, ii, iii, avr, avl, avf, v1..v6}_mv (12 leads) | ventricular |
| ptfv1_uv_ms | atrial |
| p_axis_deg | atrial |

That makes 23 secondary tests. Multiplicity is controlled with Benjamini-Hochberg
at a 5% FDR across these 23 p values (from the group coefficient) only. The
primary endpoint is not in the BH family. Both raw p and BH q are reported.

The following registered features are not secondary endpoints, for reasons that do
not depend on the data:

- `sdnn_ms`, `pnn50_pct`, `lf_ms2`, `hf_ms2`, `lf_hf_ratio`, `sampen` and
  `dfa_alpha1` need at least 300 s of signal. PTB-XL records are 10 s, so these
  features are NaN by declaration for every record.
- `rr_intervals_ms` is list-valued. It is summarised by `mean_hr_bpm` and `rmssd_ms`.
- `iab_class` (0/1/2) is categorical. Its registry declaration says it is
  secondary and descriptive only, never an endpoint. See the binary calls below.
- `mvp_score` (0–6) is a composite whose morphology term is `iab_class`, so it
  inherits that caveat (registry docstring). It is reported descriptively only:
  median (IQR) and the distribution of risk groups, with no test.
- `fqrs_present` is binary, not continuous. It is reported as a descriptive
  percentage per group, with no test, using the records where it is defined as the
  denominator. It is NaN when QRS ≥ 120 ms, because the Das narrow-QRS criteria do
  not apply there.

## Binary threshold calls: descriptive only, no test

These are reported as n (%) per group, with no p value, CI or model:

- `iab_class`: none (0), partial (1, P ≥ 120 ms) and advanced (2), plus any IAB (≥ 1)
- PTFV1 abnormal: `ptfv1_uv_ms` ≤ −4000 µV·ms
- MVP high risk: `mvp_score` ≥ 5. The low (0–2) and intermediate (3–4) groups are
  shown as well.
- `fqrs_present` = 1, as described above.

Reason: the multilead durations carry a calibration offset against method-matched
reference ranges (`gate3_sanity.md`), and on LUDB the 120 ms call agreed with the
reference at only κ 0.24–0.26. A constant offset does not bias a continuous
between-group difference, but it does change the rate of every threshold call.

## Circularity statement

PTB-XL LAO/LAE labels were assigned by cardiologists partly from ECG P-wave
criteria (P duration, PTFV1). A case-control difference in those features is
therefore evidence of **concurrent validity** of the automated measurement against
expert labels. It is not an independent clinical finding. Features outside the
LAE criteria are not circular: QRS, QT, QTc, Tpe, Tpe/QT, HR, RMSSD, fQRS and ST.

Each feature in the results table is labelled as follows:

- **circular**: `p_duration_ms`, `ptfv1_uv_ms`, plus `iab_class` and `mvp_score`,
  which are built on P duration
- **not circular**: `n_beats`, `mean_hr_bpm`, `rmssd_ms`, `qrs_duration_ms`,
  `qt_ms`, `qtc_bazett_ms`, `qtc_fridericia_ms`, `tpe_ms`, `tpe_qt_ratio`, all
  `st_j60_*`, `fqrs_present`
- **P-wave-derived, not an LAE criterion**: `p_axis_deg`. It is not one of the
  standard LAE criteria, but it is measured on the same P waves the cardiologist
  read. It is labelled separately and is not claimed as independent.

## Sensitivity analyses: planned here, not run in this step

1. **Complete matched sets** (159 sets of 1 case + 2 controls). For each set, take
   d = case value − mean of the two control values, and report the mean d with a
   95% CI (paired t), plus a Wilcoxon signed-rank p. The endpoints are the same:
   primary uncorrected, secondaries BH within their own family. Age and sex are
   matched within each set, so no covariates are used.
2. **Design B** (patients with AF in another record vs. never-AF controls). The
   records go through the identical extraction and QC (`qc.py`, current
   `config.yaml`), then the identical primary and secondary models, effect sizes
   and BH procedure. Design B labels are not P-wave-based, so the circularity
   caveat does not apply to them.

## Not done

No predictive model, classifier or decision tree is fitted. Nothing in the
delineator, feature code, `qc.py` or `config.yaml` changes for this analysis.

## Deviations

(None at the time of writing. Any later deviation is added below with its date and
reason. The text above is not edited.)

### 2026-10-10: POST-HOC check, report provenance (not pre-specified)

Added after the primary analysis (`group_comparison.md`) had been run, so it is
post-hoc and labelled as such everywhere it is reported. Reason: PTB-XL reports
were partly autogenerated by the recording device and only partly validated by
cardiologists. From `ptbxl_database.csv`, the counts of `validated_by_human`,
`initial_autogenerated_report` and `second_opinion` are reported for Design A
cases and controls separately. The primary model (`p_duration_ms ~ group + age +
sex`, OLS) is repeated on the QC-passed records whose report was human-validated,
provided that subset has at least 50 cases; otherwise only the count is reported.
Script: `paper/notes/sensitivity.py`. Results: `paper/notes/sensitivity.md`.

Sensitivity analyses 1 and 2 were run as planned above. One detail the plan does not
fix: in Sensitivity 1, BH is applied separately to the 23 paired-t p values and to
the 23 Wilcoxon p values.
