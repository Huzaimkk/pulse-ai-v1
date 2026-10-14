# Record-level QC exclusions, Design A (cohort flow)

Run 2026-09-30 with the current `config.yaml`: `pulse_ecg/qc.py`, driver
`paper/notes/qc_exclusions.py`. Per-record output `out/qc_designA.csv` (all 981
records, with `qc_pass`, every failed rule and its detail); flow table
`out/qc_flow_designA.csv`.

**QC is applied blind to case/control status.** `qc.evaluate()` receives only the
record's feature row and median-beat fiducials; the group label is not an input
(enforced by `tests/test_qc.py::test_evaluate_is_blind_to_group`), and it is joined
back only after every pass/fail decision has been made.

## Rules

| rule | fails when | source |
|---|---|---|
| delineation | the delineator reports failure | `extract()` delineation_reason |
| missing:\<feature\> | a required feature is NaN (reason carried through) | required = the four range-checked features + primary endpoint |
| missing:p_offset / qrs_onset | fiducial absent on the median beat | |
| qrs_onset_before_p_offset | qrs_onset < p_offset on the median beat | physiology; `qrs_onset_overlap.md` |
| range:hr | mean HR outside 40–200 bpm | config `qc:` |
| range:qrs_ms | QRS duration outside 50–200 ms | config `qc:` |
| range:qt_ms | QT outside 250–600 ms | config `qc:` |
| range:p_ms | P duration outside 50–200 ms | config `qc:` |

The config ranges are checked against the registry at build time. Each one must
contain its feature's declared normal range and abnormal threshold, so QC cannot
remove the finding under study (e.g. P ≥ 120 ms, QRS ≥ 120 ms). Bounds are inclusive.
The overlap rule allows qrs_onset = p_offset (PR segment 0).

## Exclusions

*Any* = records failing the rule; a record can fail several rules, so these don't sum.
*Flow* = each excluded record is attributed to its first failed rule in the order
above; these sum to the total excluded and form the flow diagram.

| rule | all: any | all: flow | case: any | case: flow | control: any | control: flow |
|---|---|---|---|---|---|---|
| Design A records | 981 | 981 | 327 | 327 | 654 | 654 |
| delineation failed | 0 | 0 | 0 | 0 | 0 | 0 |
| required feature / fiducial NaN | 0 | 0 | 0 | 0 | 0 | 0 |
| qrs_onset_before_p_offset | 152 (15.5 %) | 152 (15.5 %) | 40 (12.2 %) | 40 (12.2 %) | 112 (17.1 %) | 112 (17.1 %) |
| range:hr | 0 | 0 | 0 | 0 | 0 | 0 |
| range:qrs_ms | 51 (5.2 %) | 33 (3.4 %) | 26 (8.0 %) | 21 (6.4 %) | 25 (3.8 %) | 12 (1.8 %) |
| range:qt_ms | 10 (1.0 %) | 4 (0.4 %) | 7 (2.1 %) | 2 (0.6 %) | 3 (0.5 %) | 2 (0.3 %) |
| range:p_ms | 26 (2.7 %) | 17 (1.7 %) | 10 (3.1 %) | 5 (1.5 %) | 16 (2.4 %) | 12 (1.8 %) |
| **excluded** | **206 (21.0 %)** | **206 (21.0 %)** | **68 (20.8 %)** | **68 (20.8 %)** | **138 (21.1 %)** | **138 (21.1 %)** |
| **final n** | **775** | | **259** | | **516** | |

Percentages are of each group's starting n. No record errored, and none was dropped
without a named rule: 206 excluded + 775 kept = 981.

## Observations (reported, not acted on)

- The overlap count reproduces the independent measurement in
  `qrs_onset_overlap.md` exactly (152: 40 cases, 112 controls).
- The overall exclusion rate is balanced (20.8 % vs 21.1 %), but its makeup is not.
  Controls lose more to the overlap rule (17.1 % vs 12.2 %). Cases lose more to the
  QRS range (8.0 % vs 3.8 % any; 6.4 % vs 1.8 % as the first failure), so records
  with QRS > 200 ms or < 50 ms are enriched in LAO/LAE cases.
- Matching after QC: of 327 matched sets, 259 keep their case, 248 keep the case and
  at least one control, and 159 keep the case and both controls. The analysis plan
  must say whether broken sets are analysed unmatched, as 1:1, or dropped.

## Cohort comparability after QC — Table I draft

QC breaks the matched sets (above), so the analysis treats the groups as unmatched
and adjusts for age and sex. This table checks that the groups are still
comparable on those two variables after QC. Driver `paper/notes/gate3_sanity.py`,
output `out/table1_postqc_designA.csv`. Sex is PTB-XL `sex == 0` (male).

**Table I. Cohort characteristics, Design A (PTB-XL, sinus rhythm)**

| | Cases (LAO/LAE) | Controls | SMD |
|---|---|---|---:|
| **After QC (analysed)**, n | 259 | 516 | |
| Age, years, mean ± SD | 68.3 ± 11.9 | 68.2 ± 12.1 | 0.006 |
| Male, n (%) | 152 (58.7) | 310 (60.1) | −0.028 |
| *Before QC (matched 1:2)*, n | 327 | 654 | |
| Age, years, mean ± SD | 68.3 ± 12.2 | 68.3 ± 12.2 | 0.000 |
| Male, n (%) | 197 (60.2) | 394 (60.2) | 0.000 |

SMD = standardised mean difference, cases − controls. For age:
(mean₁ − mean₂) / √((SD₁² + SD₂²)/2). For male sex: (p₁ − p₂) / √((p₁(1−p₁) + p₂(1−p₂))/2).

Both post-QC SMDs are well below the conventional 0.1 imbalance threshold
(|0.006| and |0.028|). QC exclusions did not unbalance age or sex: the mean age
difference is 0.07 years and the male share differs by 1.4 percentage points.
Before QC the SMDs are exactly 0, as expected from exact sex matching and age
matching. Only age and sex are covered here, because they are the only covariates
in `out/cohort.csv` and the only planned adjustment variables.
