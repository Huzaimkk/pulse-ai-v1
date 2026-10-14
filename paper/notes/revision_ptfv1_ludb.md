# POST-HOC revision analysis: PTFV1 against LUDB manual V1 annotations

**This analysis is post hoc. It was added at revision and is not pre-specified in `analysis_plan.md`.** `pulse_ecg/` was not modified.

- Script: `paper/notes/revision_ptfv1_ludb.py`
- Summary: `out/revision_ptfv1_ludb.csv`
- Per-record detail: `out/revision_ptfv1_ludb_records.csv`
- Per-beat detail: `out/revision_ptfv1_ludb_beats.csv`
- PTB-XL control components: `out/revision_ptfv1_ptbxl_controls.csv`

This closes the gap noted in `gate3_sanity.md`: "The V1 single-lead boundary error has not been quantified."

## Method

- **Our V1 boundaries.** `wavedet.delineate_lead` runs on the preprocessed V1 lead with the pipeline's R-peaks. For every record, the script asserts that the per-record median equals `Delineation.per_lead["V1"]`, which is exactly what `ptfv1_uv_ms` reads.
- **LUDB V1 boundaries.** These are the cardiologist's P onset and P offset in lead V1. As in Gate 2, a P wave is paired with the detected QRS it precedes by ≤ 400 ms. Only beats inside LUDB's annotated span are scored.
- **Per-record positions.** Each position is the median, over beats, of (mark − R-peak). This is the same construction the median beat uses.
- **PTFV1, two ways, on the same median beat** (preprocessed, re-referenced to the pre-QRS baseline), using the same `atrial.terminal_negative()`:
  - (a) the registered feature on our V1 boundaries. Recomputing it reproduces the feature value for every record (asserted).
  - (b) LUDB's V1 boundaries, with nothing else changed.
- **Threshold.** The config defines abnormal as ≤ −4000 µV·ms. The request said < −4000. No value in LUDB or in the PTB-XL controls equals −4000 exactly, so the two definitions give identical results.
- **Statistics.**
  - Bland–Altman bias = mean(ours − LUDB); LoA = bias ± 1.96 SD.
  - ICC(2,1) = two-way random, absolute agreement, single measure.
  - Per-beat statistics treat beats as independent, but beats are clustered within records.

**Subsets**

- **Gate 2 held-out half, sinus** (`gate2.split_seed` 42, `tune_frac` 0.5, as in `scripts/02_validate_delineation.py`): 90 records. These are 2, 7, 11, 12, 15, 16, 17, 18, 19, 20, 22, 23, 24, 31, 32, 36, 37, 41, 42, 43, 45, 46, 48, 50, 54, 58, 59, 60, 61, 62, 64, 65, 66, 70, 72, 73, 74, 75, 77, 78, 79, 84, 91, 92, 94, 98, 100, 106, 111, 114, 117, 118, 120, 121, 123, 126, 130, 137, 138, 140, 141, 142, 143, 146, 150, 152, 155, 157, 159, 161, 162, 163, 164, 167, 168, 169, 170, 176, 178, 180, 182, 183, 184, 185, 186, 193, 194, 195, 196, 199.
- **All sinus:** 182 records.
- **All 200 records.** Each of the 18 non-sinus records lacks either our V1 P marks or LUDB's, so none of them enters the paired comparisons. As a result, the all-200 paired results are identical to the all-sinus ones; only the detection and % abnormal rows differ.

## 1. V1 boundary agreement (ours − LUDB, ms)

| Subset | Boundary | Level | n | Bias | SD | LoA | MAE |
|---|---|---|---|---|---|---|---|
| Held-out sinus | P onset | record | 86 | −2.3 | 42.2 | −85.1 to 80.5 | 29.8 |
| | P onset | beat | 632 | −7.2 | 43.8 | −93.0 to 78.6 | 32.3 |
| | P offset | record | 86 | +37.1 | 40.4 | −42.1 to 116.2 | 42.4 |
| | P offset | beat | 632 | +34.2 | 39.7 | −43.5 to 111.9 | 40.7 |
| All sinus / all 200 | P onset | record | 173 | −0.7 | 41.9 | −82.9 to 81.4 | 30.1 |
| | P onset | beat | 1214 | −3.7 | 43.5 | −88.9 to 81.6 | 31.7 |
| | P offset | record | 173 | +42.0 | 37.7 | −31.8 to 115.7 | 46.8 |
| | P offset | beat | 1214 | +39.4 | 38.4 | −35.9 to 114.8 | 45.1 |

**V1 P detection**

| Subset | Level | Ours detects when LUDB has V1 P | Ours detects when LUDB has none |
|---|---|---|---|
| Held-out sinus | record | 86/87 (98.9 %) | 3/3 (100 %) |
| | beat | 632/707 (89.4 %) | 88/110 (80.0 %) |
| All sinus | record | 173/176 (98.3 %) | 6/6 (100 %) |
| | beat | 1214/1398 (86.8 %) | 171/223 (76.7 %) |
| All 200 | record | 173/176 (98.3 %) | 23/24 (95.8 %) |
| | beat | 1214/1398 (86.8 %) | 300/432 (69.4 %) |

Our method places a P wave in V1 in almost every record where LUDB annotates none. Most of those records are non-sinus (all 200: 23 of 24), so they are outside the pipeline's population, but this is still a false-positive tendency.

## 2. PTFV1 agreement: (a) our V1 boundaries vs (b) LUDB V1 boundaries, same signal

| Subset | n | Bias (a − b), µV·ms | SD | LoA | Spearman r | ICC(2,1) | κ (≤ −4000) | Abnormal: a / b / both / a only / b only |
|---|---|---|---|---|---|---|---|---|
| Held-out sinus | 86 | −605 | 1218 | −2993 to 1782 | 0.925 | 0.928 | 0.785 | 25 / 18 / 18 / 7 / 0 |
| All sinus (= all 200) | 173 | −725 | 1264 | −3203 to 1753 | 0.933 | 0.896 | 0.774 | 47 / 33 / 33 / 14 / 0 |

**Components (a − b), all sinus, n = 173**

| Component | Bias | SD | LoA |
|---|---|---|---|
| Terminal-negative duration, ms | +17.2 | 28.0 | −37.6 to 72.0 |
| Terminal-negative depth, µV | +0.8 | 5.8 | −10.5 to 12.2 |
| Held-out terminal-negative duration, ms | +14.3 | 27.8 | −40.2 to 68.7 |

The disagreement is one-directional. Our boundaries never call a record normal that LUDB's boundaries call abnormal (0 such records). They make PTFV1 more negative because the late V1 P offset (+37 to +42 ms) lengthens the terminal-negative run. The terminal-negative depth is essentially unchanged.

## 3. % abnormal PTFV1 (≤ −4000 µV·ms)

| Population | Boundaries | n | Abnormal | Median (IQR), µV·ms |
|---|---|---|---|---|
| LUDB held-out sinus | ours | 89 | 25 (28.1 %) | −1629 (−4247 to 0) |
| | LUDB | 87 | 18 (20.7 %) | −594 (−3096 to 0) |
| LUDB all sinus | ours | 179 | 47 (26.3 %) | −1629 (−4235 to 0) |
| | LUDB | 176 | 33 (18.8 %) | −534 (−3009 to 0) |
| LUDB all 200 | ours | 196 | 49 (25.0 %) | −1504 (−3895 to 0) |
| | LUDB | 176 | 33 (18.8 %) | −534 (−3009 to 0) |
| **PTB-XL analysed controls** (`group_threshold_calls_designA.csv`) | ours | 516 | **258 (50.0 %)** | −4003 (−5925 to −1975) |

LUDB is a hospital population, not a healthy reference, so it gives context, not a normal rate. On LUDB, the same pipeline calls 26 % abnormal, against 19 % with the manual boundaries. The PTB-XL controls reach 50 %, about twice the LUDB rate from the same code.

## 4. Is the 50 % control rate a sign, units or threshold error?

| Check | Result |
|---|---|
| Recomputed vs stored control PTFV1 | 516/516 identical |
| Values > 0 (impossible by definition) | 0 |
| Values = 0 (upright V1 P, no terminal negative component) | 60 |
| Values exactly −4000 (≤ vs < ambiguity) | 0, so % ≤ −4000 = % < −4000 = 50.0 % |
| Control PTFV1 median (IQR) | −4003 (−5925 to −1975) µV·ms |
| Units | Signals are in mV (`io._check_units` rejects any other unit). `terminal_negative` returns depth = −min × 1000 in µV and duration = samples × 1000/fs in ms, so the product is µV·ms. The threshold −4000 µV·ms = 0.04 mm·s at 10 mm/mV (Morris) |
| Control terminal-negative duration, median (IQR) | 74 (60–94) ms |
| Control terminal-negative depth, median (IQR) | 50.0 (31.9–68.2) µV |
| Control V1 single-lead P duration, median (IQR) | 133 (109–152) ms |
| LUDB sinus, for comparison (medians) | terminal duration 58 ms on our boundaries vs 18 ms on LUDB's; depth 26.3 vs 26.7 µV |

There is no sign error: no control has a positive value, and negative values mean a terminal negative component, as Morris defines it. There is no units error: the inputs are checked to be in mV and the formula yields µV·ms. There is no threshold error: −4000 µV·ms is the Morris 0.04 mm·s cut, and the ≤/< choice moves no record. The 50 % falls out of the distribution itself: the control median, −4003, sits almost exactly on the threshold.

Two measurable contributors push the controls' distribution to that point:
- **Terminal-negative duration.** The controls' median (74 ms) is long, and section 2 shows our V1 offset adds about 14–17 ms of terminal duration relative to manual boundaries.
- **Terminal-negative depth.** Median depth in PTB-XL controls is about 50 µV, against about 26 µV in LUDB. The boundaries do not affect depth, so this difference comes from the PTB-XL signals or population, not from boundary placement.

## Verdict (post hoc)

**Limited.** As a continuous measurement, PTFV1 ranks records in close agreement with PTFV1 computed from cardiologist V1 boundaries on the same signal: Spearman 0.93, ICC(2,1) 0.93 on the held-out sinus half and 0.90 on all sinus records. This supports the pre-specified use of PTFV1 as a continuous secondary outcome. The V1 boundaries themselves are not validated. P offset is late by +37 to +42 ms, and onset and offset limits of agreement are about ±80 ms, far wider than inter-observer tolerances for P boundaries. That error is systematic. It makes PTFV1 more negative by about 600–700 µV·ms on average and turns 7 of 25 held-out abnormal calls into calls the manual boundaries do not make (κ 0.78). No disagreement goes the other way. The binary "abnormal PTFV1" call therefore overestimates prevalence and should stay descriptive only, as the analysis plan already specifies.

The 50 % abnormal rate in PTB-XL controls is not a sign, units or threshold error. It reflects a control distribution centred on the threshold, from long terminal-negative durations, which are partly inflated by the late V1 offset, and from terminal-negative depths about twice those in LUDB. This analysis cannot say whether those depths reflect the population or the recording. The rate should not be presented as a population prevalence.

The case–control difference in PTFV1 was computed with the same boundary method in both groups. It is not invalidated by a common bias, but this analysis does not show that the bias is equal in cases and controls.
