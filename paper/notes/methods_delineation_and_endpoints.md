# Methods notes — delineation (Gate 2) and endpoints

Decisions agreed 29 Sep 2026. Numbers come from `scripts/02_validate_delineation.py`
(`out/gate2_summary.csv`, `out/gate2_p_duration.csv`) and `scripts/01_cohort.py`
(`out/gate1.csv`, `out/table1.csv`). Re-run those scripts before any number enters the paper (Rule 8).

## Cohort (Gate 1, Table I)

PTB-XL v1.0.3, records500. Every analysed record is in sinus rhythm; one record per
patient; controls matched exactly on sex and on age within ±5 years, 2 per case (seed 42).
SCP likelihood 0 means "stated", so it counts; PTB-XL rhythm statements are always
likelihood 0, so the likelihood sensitivity threshold applies to diagnostic labels only.

- **Design A (primary):** `LAO/LAE` cases vs matched controls — 327 cases, 654 controls,
  mean age 68.3, 60.2 % male in both groups. Gate 1 passes.
- **Design B (sensitivity analysis):** sinus records of patients with AFIB in another
  record — 118 cases, 236 controls, mean age 70.2.
- At likelihood 100, design A keeps 230 cases; design B collapses to 5 cases because
  1,466 of 1,514 AFIB statements carry likelihood 0, so that check is uninformative for B.

## Delineation

Wavelet delineator of Martínez et al. [martinez2004wavedet], each lead delineated
separately, then combined per beat with the k-of-N agreement rule of Laguna et al.
[laguna1994multilead]. Code: `pulse_ecg/wavedet.py`; parameters: `config.yaml` → `wavedet:`.

State these deviations from the published method:
1. QRS detection uses our multilead R-peak detector, not the paper's modulus-maximum-line detector.
2. At 500 Hz the same quadratic-spline filters are used with the scale index shifted by one
   (same analogue bands); the paper resamples its 250 Hz filter bank instead.
3. Search windows are ours (the paper does not specify them); significant P/T slopes must lie
   within 100 / 150 ms of the dominant slope.
4. A single-lead P wave needs a deflection ≥ 0.02 mV, and a global mark needs k + 1 leads.
   Without this, ε_P (A.4) accepted QRS/T leakage as a P wave when none existed.
5. Multilead rule: k = 3 with δ from Laguna 1994 for QRS and T; for P boundaries k = 2, δ = 10 ms.
   With k = 3 / δ = 6 ms, P duration came out 17 ms short on the tune half.

neurokit2's DWT delineator was evaluated first and rejected: on the tune half, QRS onset was
−115 ms and P onset −31 ± 53 ms (mean ± SD).

### Validation protocol

- LUDB split 50/50 by record (seed 42). Every tunable parameter was tuned on the tune half
  only; the figures reported are from the held-out half, sinus-rhythm records (the population
  analysed). LUDB annotates only the middle of each record; only reference beats inside the
  annotated span are scored.
- QTDB is a fully external check: 2 leads, 250 Hz (4 ms per sample vs 2 ms for LUDB), errors in ms.
- The within-record SD (the mean of per-record SDs) is the statistic of Martínez 2004,
  reported beside the pooled SD. Martínez chose the best channel per point against the
  cardiologist on QTDB, so their figures are optimistic relative to any real multilead rule.
- CSE 2-SD tolerances [cse1985standards], as tabulated in Martínez 2004 Table III:
  P on 10.2, P end 12.7, QRS on 6.5, QRS end 11.6, T end 30.6 ms.

### The LUDB reference and P duration

LUDB annotates each lead independently. The earliest onset / latest offset over 12 leads
is an order statistic that inflates the span by construction: 80 % of LUDB sinus records
then have P ≥ 120 ms, against a plausible IAB prevalence of 10–20 %. We made one attempt
at a fairer reference, pre-specified as the same k-of-N rule used for measurement, and
report both constructions without further tuning.

Held-out LUDB, sinus records, P duration (per-record median):

| Reference | Record-level error, mean ± SD (ms) | MAE (ms) | Ref. prevalence ≥ 120 ms | κ at 120 ms |
|---|---|---|---|---|
| min/max over leads | −8.0 ± 28.8 | 20.4 | 80 % | 0.24 |
| k-of-N rule | +5.7 ± 23.4 | 16.3 | 66 % | 0.26 |
| QTDB (external, 1 annotator) | +3.3 ± 19.6 | 15.1 | 40 % | 0.34 |

κ did not improve materially. State the problem plainly: agreement on the 120 ms dichotomy
is low. It reflects both our measurement error and the reference itself; even the
agreement-rule reference puts two thirds of LUDB sinus records above 120 ms. P duration is
therefore not used as a per-record classifier.

## Endpoints

- **Primary:** P-wave duration as a continuous measure — design A cases vs matched controls.
  Report the mean difference in ms and a standardised effect size, each with a 95 % CI.
  Dichotomising a continuous variable at a threshold loses power and is widely criticised.
- **Secondary, descriptive:** IAB classification (P ≥ 120 ms, morphology), reported with the
  measured delineation error above.
- Design B repeats the primary analysis as a sensitivity analysis.

PTFV1 and P axis are unaffected by the multilead pooling problem. PTFV1 is measured in V1
alone, with no multilead pooling; its duration term still uses V1's single-lead P boundaries.
P axis is an area ratio (I vs aVF), which is less sensitive to exact boundaries. The atrial
axis stands on these two.

Anticoagulation: never recommend or imply it (CLAUDE_1.md Rule 9; ARCADIA, [kamel2024arcadia]).

## Known limitation of the synthetic fixture

Six synthetic delineation tests (`tests/test_delineate.py::test_boundaries_match_truth`)
fail: wavedet places P onset 8–18 ms early and QRS/T boundaries about 10 ms outward on the
synthetic ECG. The fixture builds each wave from half-sines with exact support, so every
wave starts and ends with a slope discontinuity that no physiological signal has; the
wavelet scales respond to that kink. This is a limitation of the fixture, not of the
delineator, whose error is measured on cardiologist-annotated data above. The tests are
deliberately left failing and must not be loosened to pass.

## Next (after /clear)

Build the features in order — rate → ventricular → atrial → quality — with tests written
before each module (Code_Architecture_Plan.md §7), then Gate 3 (every index in its
published normal range, via `REGISTRY.sanity`). The `extract()` entry point and
`report.prf` that `tests/test_features.py` imports do not exist yet.
