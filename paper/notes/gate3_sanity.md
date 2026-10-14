# Gate 3 — sanity check, Design A after QC

Run on the **775 QC-passed Design A records** (259 cases, 516 controls;
`qc_exclusions.md`). Driver `paper/notes/gate3_sanity.py`, which reads
`out/qc_designA.csv` and `out/cohort.csv` only, so nothing is re-extracted.
Outputs: `out/gate3_sanity_designA.csv` (all features) and
`out/gate3_offsets_designA.csv` (residual offsets). Values are pooled over cases
and controls. No outcome feature is split by group here.

Revised 2026-10-10. The first version (2026-09-30) judged P and QRS duration
against single-lead visual-caliper ranges. The comparison is corrected below. The
extraction is unchanged.

## Why the reference ranges changed

The classic ranges, P 80–120 ms and QRS 70–110 ms, come from single-lead caliper
reading. AHA/ACCF/HRS Part I [kligfield2007aha1] states that single-lead
measurement systematically underestimates wave durations, that measurement on
simultaneous leads is systematically greater, and that population criteria for P
and QRS duration need redefinition for it. PULSE measures both durations
multilead: global onset to global offset across all simultaneous leads
(`p_duration_ms` on the 12-lead set; `methods_delineation_and_endpoints.md`). So
the correct comparator is a lead-independent, age-matched healthy population.

**Comparator.** Rijnbeek et al. [rijnbeek2014normal] give lead-independent normal
values for 13,354 healthy subjects. Their medians at age 60–89 are P 120–121 ms
(men) and 114–118 ms (women), and QRS 98–101 ms (men) and 92 ms (women). For this
cohort (mean age 68, 60 % male) that gives **P ≈ 115–125 ms** and
**QRS ≈ 95–100 ms**. The age dependence of P-wave duration in healthy subjects
[havmoller2007pwave] and in a large unselected population [giovanardi2022ecg],
and P-duration distributions in a large primary-care cohort [nielsen2015pwave],
support using an age-matched comparator rather than the young-adult caliper range.

## Gate 3 indices

| feature | n | median | IQR | comparator | residual offset (median − comparator) | outcome |
|---|---:|---:|---|---|---|---|
| P duration (ms), primary endpoint | 775 | 143 | 133–156 | 115–125 (method + age matched) | **+18 to +28 ms** | calibration offset |
| QRS duration (ms) | 775 | 112 | 100–130 | 95–100 (method + age matched) | **+12 to +17 ms** | calibration offset |
| PR interval (ms) † | 775 | 186 | 170–204 | 120–200 | in range (29.4 % > 200) | PASS |
| QT (ms) | 775 | 426 | 396–468 | 350–450 | in range (34.2 % > 450) | PASS |
| mean HR (bpm) | 775 | 72.1 | 64.6–82.6 | 60–100 | in range (13.3 % outside) | PASS |
| Tpe/QT | 766 | 0.265 | 0.235–0.327 | none established | — | runs high; **secondary** |

† PR is not a registered feature. It is derived from the same median-beat
fiducials as `p_duration_ms`: PR = qrs_onset − p_onset. It therefore carries the
same P-onset error as P duration.

**Sex-matched check (age 60–89, pooled over cases and controls).** The offset is
about the same in men and women. That is what a calibration shift would produce,
rather than a problem confined to one subgroup.

| feature | stratum | n | median | Rijnbeek median | residual offset |
|---|---|---:|---:|---|---|
| P duration (ms) | men 60–89 | 303 | 146 | 120–121 | +25 to +26 ms |
| P duration (ms) | women 60–89 | 245 | 142 | 114–118 | +24 to +28 ms |
| QRS duration (ms) | men 60–89 | 303 | 116 | 98–101 | +15 to +18 ms |
| QRS duration (ms) | women 60–89 | 245 | 107 | 92 | +15 ms |

### What the offsets mean

- **P duration (+18 to +28 ms).** A constant offset is a calibration bias, not a discrimination failure. It shifts cases and controls equally, so the continuous group comparison stands, while absolute values and threshold-based calls do not.
- **QRS duration (+12 to +17 ms).** A constant offset is a calibration bias, not a discrimination failure. It shifts cases and controls equally, so the continuous group comparison stands, while absolute values and threshold-based calls do not.

Limits on that argument, stated so the paper does not overclaim:

- **"Constant" and "equally" are assumptions, not results of this note.** Testing them needs a group comparison, which has not been run.
- **The P-duration error also has a slope.** On synthetic ECGs, measured = 0.782 × true + 27.9 ms (`p_duration_compression.md`). The slope shrinks a true case-control difference to 0.78 of its size, so it biases the continuous comparison toward the null (conservative), not away from it.
- **There is a random per-record component too.** On LUDB the P-duration bias is +20 to +34 ms with SD about 35 ms (`gate2_ludb.md`). That costs power but adds no bias.
- **No slope has been measured for QRS.** Its offset is treated as constant by assumption.

## Threshold-based features affected by P-wave timing error

A positive bias in P-wave timing inflates every call that applies a fixed
threshold to a duration, or to a product that contains one.

| feature | threshold | why it is inflated | share called abnormal (775) |
|---|---|---|---:|
| `iab_class` | P duration ≥ 120 ms defines IAB [bayesdeluna2012iab] | the offset pushes typical P waves over 120 ms | 93.5 % (725) with IAB ≥ 1 |
| `ptfv1_uv_ms`, abnormal call | ≤ −4000 µV·ms [morris1964ptfv1] | PTFV1 = terminal duration × amplitude, so any lengthening of the measured terminal duration scales the product. Its duration comes from V1's own single-lead P boundaries (`atrial.py`), not the multilead P duration, so the +18 to +28 ms global offset does not carry over directly. The V1 single-lead boundary error has not been quantified. | 60.6 % ≤ −4000 |
| `mvp_score` | high risk ≥ 5 [alexander2019mvp] | its P-duration and IAB components are thresholded on the inflated duration | median 4 (IQR 4–5); 93.5 % above 0–2 |

**These are analysed as continuous values; their binary calls are descriptive
only.** PTFV1 enters the analysis in µV·ms. IAB and the MVP duration component
rest on P duration, which enters as the continuous primary endpoint. The
threshold calls (IAB present, PTFV1 abnormal, MVP high risk) are reported for
description and are not used as analysis endpoints.

## Tpe/QT — secondary

No age-matched, method-matched comparator was established for Tpe/QT. The
0.15–0.25 range is kept for reference only. Against it, the cohort runs high:
median 0.265 (IQR 0.235–0.327), with 60.3 % of records above 0.25. Tpe has a
median of 109 ms (IQR 96–150). The cause is the T offset, the least accurate
boundary in Gate 2: t_offset SD 40–109 ms on LUDB test against a CSE 2SD
tolerance of 30.6 ms (`out/gate2_summary.csv`). **Tpe/QT is marked secondary.**
It is not a sanity anchor for the pipeline and not a primary outcome.

## Other registered features with a declared range (unchanged)

These were not part of the correction. Same median rule against the declared range.

| feature | n | median | IQR | declared range | % outside | outcome |
|---|---:|---:|---|---|---:|---|
| QTc Bazett (ms) | 775 | 464 | 434–515 | 390–450 | 63.2 | outside: inherits the long QT tail; Bazett over-corrects at HR > 60 |
| QTc Fridericia (ms) | 775 | 449 | 422–496 | 390–450 | 53.4 | PASS (by 1 ms) |
| P axis (deg) | 775 | 55 | 40–67 | 0–75 | 12.9 | PASS |

PTFV1 and MVP score are covered in the threshold section above.

## Registered features with no declared range (reported, not gated)

| feature | n | median | IQR | note |
|---|---:|---:|---|---|
| n_beats | 775 | 12 | 11–14 | |
| rmssd_ms | 775 | 17.3 | 9.8–33.7 | |
| sdnn_ms, pnn50_pct, lf_ms2, hf_ms2, lf_hf_ratio, sampen, dfa_alpha1 | 0 | — | — | NaN by design: need ≥ 300 s, PTB-XL is 10 s |
| tpe_ms | 766 | 109.3 | 96.0–150.0 | 9 NaN |
| fqrs_present | 485 | 0 | 0–0 | 67/485 positive; NaN where QRS ≥ 120 ms (Das narrow-QRS criteria) |
| st_j60_i_mv | 775 | −0.008 | −0.039 to 0.018 | |
| st_j60_ii_mv | 775 | 0.012 | −0.018 to 0.038 | |
| st_j60_iii_mv | 775 | 0.017 | −0.004 to 0.043 | |
| st_j60_avr_mv | 775 | −0.000 | −0.025 to 0.024 | |
| st_j60_avl_mv | 775 | −0.012 | −0.038 to 0.008 | |
| st_j60_avf_mv | 775 | 0.014 | −0.007 to 0.036 | |
| st_j60_v1_mv | 775 | 0.061 | 0.028 to 0.102 | |
| st_j60_v2_mv | 775 | 0.111 | 0.059 to 0.180 | |
| st_j60_v3_mv | 775 | 0.090 | 0.041 to 0.160 | |
| st_j60_v4_mv | 775 | 0.037 | −0.004 to 0.087 | |
| st_j60_v5_mv | 775 | 0.002 | −0.042 to 0.039 | |
| st_j60_v6_mv | 775 | −0.006 | −0.047 to 0.022 | |
| iab_class | 775 | 1 | 1–2 | 0: 50, 1: 522, 2: 203 |

`rr_intervals_ms` is list-valued and is summarised by `mean_hr_bpm`.

## Registry notes (reported, not changed)

- The registry still declares P 80–120 ms and QRS 70–110 ms, the single-lead ranges. They feed `REGISTRY.sanity()` and Table III. As instructed, the registry and `config.yaml` are unchanged. The method-matched comparators live in `gate3_sanity.py` only.
- PR is not a registered feature, and `mean_hr_bpm` declares no normal range.

## Gate 3 outcome

- **P and QRS duration** fall outside the single-lead ranges but carry a residual positive calibration offset against method- and age-matched norms: +18 to +28 ms (P) and +12 to +17 ms (QRS), about equal in men and women.
- **PR, QT and mean HR** are in range.
- **Tpe/QT** runs high and is secondary.

On this basis the continuous case-control comparison of P duration can proceed. It
has not been started. Absolute durations and threshold-based calls are not
interpreted as calibrated values.
