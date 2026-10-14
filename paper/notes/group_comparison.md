# Case-control comparison, Design A (Step B)

Run 2026-10-10. The analysis follows the pre-specified plan in
`paper/notes/analysis_plan.md` (sha256 `fb96c390…ad1b0`, recorded by the script at
run time), and that file has not been edited since. Script:
`paper/notes/group_comparison.py`. Outputs: `out/group_comparison_designA.csv`
and `out/group_threshold_calls_designA.csv`. Records: 775 QC-passed Design A
records (259 cases, 516 controls).
No deviations from the plan.

## Primary result

**P-wave duration, adjusted for age and sex (cases minus controls): +2.0 ms
(95% CI −0.9 to 5.0), p = 0.17.** Hedges' g = 0.10 (−0.05 to 0.25), Cliff's
δ = 0.10.

The pre-specified primary test is not significant. The CI rules out an adjusted
difference larger than about 5 ms on the measured scale.

## All features

The model for each feature is `feature ~ group + age + sex` (OLS). Diff is the
adjusted difference, cases minus controls. Hedges' g and Cliff's δ are unadjusted.
q is the Benjamini-Hochberg value across the 23 secondary endpoints only; the
primary endpoint is not in that family. The circularity labels follow the plan:
*circular* means the feature is part of the ECG LAE criteria that PTB-XL
cardiologists used to assign the case label. ST values are in mV and PTFV1 is in
µV·ms.

| feature | circularity | n case / control | case median (IQR) | control median (IQR) | adj. diff (95% CI) | p | Hedges' g (95% CI) | Cliff's δ | BH q |
|---|---|---|---|---|---|---|---|---|---|
| **p_duration_ms** (primary) | circular | 259 / 516 | 146.0 (136.0–156.0) | 141.0 (131.0–156.0) | 2.04 (−0.88 to 4.95) | 0.1701 | 0.10 (−0.05 to 0.25) | 0.10 | — |
| ptfv1_uv_ms | circular | 258 / 516 | −7607 (−10757 to −5090) | −4003 (−5925 to −1975) | −3627 (−4237 to −3017) | <0.0001 | −0.88 (−1.04 to −0.73) | −0.47 | <0.0001 |
| p_axis_deg | P-derived, not LAE criterion | 259 / 516 | 57.2 (45.8–67.9) | 54.5 (36.3–66.4) | 4.87 (0.34 to 9.40) | 0.0352 | 0.16 (0.01 to 0.31) | 0.12 | 0.0405 |
| n_beats | not circular | 259 / 516 | 13 (11–14) | 12 (11–13) | 0.84 (0.54 to 1.14) | <0.0001 | 0.42 (0.27 to 0.57) | 0.23 | <0.0001 |
| mean_hr_bpm | not circular | 259 / 516 | 76.6 (67.5–86.8) | 70.2 (64.0–79.2) | 4.94 (3.15 to 6.74) | <0.0001 | 0.41 (0.26 to 0.56) | 0.23 | <0.0001 |
| rmssd_ms | not circular | 259 / 516 | 16.2 (8.3–37.4) | 17.6 (10.4–32.7) | 2.63 (−9.93 to 15.2) | 0.6817 | 0.03 (−0.12 to 0.18) | −0.04 | 0.7127 |
| qrs_duration_ms | not circular | 259 / 516 | 120.0 (107.0–148.0) | 108.0 (99.0–124.0) | 14.2 (10.4 to 17.9) | <0.0001 | 0.54 (0.39 to 0.70) | 0.30 | <0.0001 |
| qt_ms | not circular | 259 / 516 | 442.0 (408.0–488.0) | 418.0 (390.8–458.0) | 20.5 (13.1 to 27.9) | <0.0001 | 0.41 (0.26 to 0.56) | 0.24 | <0.0001 |
| qtc_bazett_ms | not circular | 259 / 516 | 487.5 (451.8–546.2) | 454.3 (425.9–497.5) | 39.6 (30.9 to 48.3) | <0.0001 | 0.67 (0.51 to 0.82) | 0.35 | <0.0001 |
| qtc_fridericia_ms | not circular | 259 / 516 | 470.3 (436.3–529.5) | 439.2 (415.4–480.1) | 32.8 (24.9 to 40.6) | <0.0001 | 0.62 (0.47 to 0.77) | 0.33 | <0.0001 |
| tpe_ms | not circular | 256 / 510 | 120.8 (100.8–174.5) | 105.8 (94.0–138.4) | 17.7 (11.7 to 23.6) | <0.0001 | 0.44 (0.29 to 0.59) | 0.24 | <0.0001 |
| tpe_qt_ratio | not circular | 256 / 510 | 0.289 (0.238–0.361) | 0.258 (0.234–0.312) | 0.026 (0.016 to 0.036) | <0.0001 | 0.37 (0.22 to 0.52) | 0.19 | <0.0001 |
| st_j60_i_mv | not circular | 259 / 516 | −0.023 (−0.066 to 0.007) | −0.002 (−0.026 to 0.023) | −0.023 (−0.031 to −0.015) | <0.0001 | −0.42 (−0.58 to −0.27) | −0.27 | <0.0001 |
| st_j60_ii_mv | not circular | 259 / 516 | 0.004 (−0.030 to 0.038) | 0.017 (−0.011 to 0.039) | −0.011 (−0.020 to −0.003) | 0.0104 | −0.20 (−0.35 to −0.05) | −0.14 | 0.0126 |
| st_j60_iii_mv | not circular | 259 / 516 | 0.024 (−0.006 to 0.055) | 0.014 (−0.004 to 0.038) | 0.012 (0.003 to 0.021) | 0.0078 | 0.20 (0.05 to 0.35) | 0.12 | 0.0100 |
| st_j60_avr_mv | not circular | 259 / 516 | 0.010 (−0.018 to 0.039) | −0.007 (−0.028 to 0.018) | 0.017 (0.010 to 0.024) | <0.0001 | 0.36 (0.21 to 0.51) | 0.24 | <0.0001 |
| st_j60_avl_mv | not circular | 259 / 516 | −0.020 (−0.058 to 0.000) | −0.007 (−0.030 to 0.011) | −0.018 (−0.025 to −0.010) | <0.0001 | −0.35 (−0.50 to −0.20) | −0.23 | <0.0001 |
| st_j60_avf_mv | not circular | 259 / 516 | 0.013 (−0.014 to 0.042) | 0.015 (−0.004 to 0.035) | 0.001 (−0.007 to 0.008) | 0.8820 | 0.01 (−0.14 to 0.16) | −0.03 | 0.8820 |
| st_j60_v1_mv | not circular | 259 / 516 | 0.070 (0.028 to 0.121) | 0.057 (0.026 to 0.091) | 0.029 (0.014 to 0.044) | 0.0002 | 0.29 (0.14 to 0.44) | 0.14 | 0.0002 |
| st_j60_v2_mv | not circular | 259 / 516 | 0.127 (0.054 to 0.206) | 0.106 (0.059 to 0.168) | 0.047 (0.024 to 0.069) | <0.0001 | 0.30 (0.15 to 0.45) | 0.08 | <0.0001 |
| st_j60_v3_mv | not circular | 259 / 516 | 0.098 (0.044 to 0.197) | 0.087 (0.041 to 0.151) | 0.041 (0.020 to 0.062) | 0.0002 | 0.28 (0.13 to 0.43) | 0.09 | 0.0002 |
| st_j60_v4_mv | not circular | 259 / 516 | 0.027 (−0.019 to 0.089) | 0.040 (0.005 to 0.086) | −0.003 (−0.018 to 0.011) | 0.6379 | −0.04 (−0.19 to 0.11) | −0.10 | 0.6987 |
| st_j60_v5_mv | not circular | 259 / 516 | −0.028 (−0.079 to 0.023) | 0.011 (−0.024 to 0.044) | −0.038 (−0.049 to −0.026) | <0.0001 | −0.48 (−0.63 to −0.33) | −0.30 | <0.0001 |
| st_j60_v6_mv | not circular | 259 / 516 | −0.034 (−0.086 to 0.008) | 0.001 (−0.028 to 0.029) | −0.039 (−0.049 to −0.029) | <0.0001 | −0.59 (−0.74 to −0.44) | −0.33 | <0.0001 |
| mvp_score (descriptive) | circular | 259 / 516 | 4 (4–5) | 4 (4–5) | not tested | — | — | — | — |

Twenty of the 23 secondary endpoints have q < 0.05. The three with q ≥ 0.05 are
RMSSD, ST aVF and ST V4.

The following registered features are not in the table, for the reasons given in
the plan: `sdnn_ms`, `pnn50_pct`, `lf_ms2`, `hf_ms2`, `lf_hf_ratio`, `sampen` and
`dfa_alpha1` are NaN on every 10 s record (the script asserts this), and
`rr_intervals_ms` is list-valued.

## Threshold calls: descriptive only, no test

| call | circularity | cases n/N (%) | controls n/N (%) |
|---|---|---|---|
| IAB none (0) | circular | 17/259 (6.6) | 33/516 (6.4) |
| IAB partial (1) | circular | 177/259 (68.3) | 345/516 (66.9) |
| IAB advanced (2) | circular | 65/259 (25.1) | 138/516 (26.7) |
| any IAB (≥ 1) | circular | 242/259 (93.4) | 483/516 (93.6) |
| PTFV1 abnormal (≤ −4000 µV·ms) | circular | 211/258 (81.8) | 258/516 (50.0) |
| MVP low (0–2) | circular | 17/259 (6.6) | 33/516 (6.4) |
| MVP intermediate (3–4) | circular | 123/259 (47.5) | 251/516 (48.6) |
| MVP high (≥ 5) | circular | 119/259 (45.9) | 232/516 (45.0) |
| fQRS present | not circular | 30/127 (23.6) | 37/358 (10.3) |

The fQRS denominators are smaller because fQRS is defined only when QRS < 120 ms
(Das criteria). More cases than controls have a wide QRS, so more case records
are NaN.

## Reading the results (interpretation; none of these points are pre-specified tests)

- **Primary: no detectable difference in P duration.** This is a null result for
  the pre-specified endpoint and should be reported as one. Two known measurement
  properties are relevant context, but they do not change the result:
  1. The delineator compresses P duration with a slope of 0.78
     (`p_duration_compression.md`). That attenuates any true difference, so on
     the true scale the CI upper bound would be about 5.0 / 0.78 ≈ 6 ms. This is
     a back-of-envelope rescaling, not part of the plan.
  2. The threshold calls show the absolute-level problem directly. 93% of both
     groups are classed as IAB, which is consistent with the calibration offset
     in `gate3_sanity.md`, and the 120 ms call does not separate the groups.
- **PTFV1 separates the groups strongly**: g −0.88, abnormal call 82% vs 50%.
  PTFV1 is a circular feature, so this is concurrent validity: the automated V1
  terminal-force measurement agrees with the cardiologists' LAE labelling. It is
  not an independent finding. A plausible reading of the contrast with the null P
  duration is that the PTB-XL LAO/LAE label relied more on V1 morphology than on
  P duration, or that global P duration is the noisier measurement. These data
  cannot tell those two explanations apart.
- **The non-circular features differ substantially.** Cases have a wider QRS
  (+14 ms), a longer QTc (Fridericia +33 ms) and Tpe, a heart rate about 5 bpm
  higher, and ST depression in the lateral leads (I, aVL, V5, V6) with ST
  elevation in aVR and V1–V3. Taken together, the pattern looks like ventricular
  structural or conduction disease travelling with the LAE label (for example
  LVH with strain). These features are independent of the label's P-wave
  criteria. That does not make them independent of the cardiologist's overall
  reading of the ECG. The model adjusts for age and sex only, and other diagnoses
  were not adjusted for.
- Bazett overcorrects at higher heart rates, and cases have a higher heart rate.
  The Fridericia difference is smaller (33 vs 40 ms) but is still large.

## Not done in this step

The sensitivity analyses (159 complete matched sets; Design B) were planned but
not run. No model, classifier or decision tree was fitted. No changes were made
to the delineator, the features, `qc.py` or `config.yaml`.
