# POST-HOC revision analysis: records included vs excluded by QC (Design A)

**This analysis is post hoc. It was added at revision and is not pre-specified in `analysis_plan.md`.** The QC rules, `out/qc_designA.csv`, `out/qc_flow_designA.csv` and every other existing output are unchanged; this analysis only reads them.

Script: `paper/notes/revision_qc_included_excluded.py`. It writes three files:
- `out/revision_qc_included_excluded.csv`: characteristics and tests
- `out/revision_qc_reasons.csv`: exclusions by rule
- `out/revision_qc_flow_chain.csv`: the cohort chain

Population: all 981 matched Design A records that entered QC (327 cases, 654 controls). 775 were included and 206 excluded.

## 1. Cohort chain (eligible → usable → matched → analysed)

| Step | Cases | Controls | Reason for loss |
|---|---|---|---|
| Eligible (sinus record, first per patient; cases carry LAO/LAE, controls never do) | 332 | 14,700 | |
| − records500 file absent | 0 | 0 | `file_present` |
| − age < 40 | 5 | 2,113 | `cohort.min_age: 40` |
| **Usable** | **327** | **12,587** | |
| − not selected by matching | 0 | 11,933 | 1:2 matching, exact on sex, age ±5 y; every usable case found 2 controls |
| **Matched (entered QC)** | **327** | **654** | |
| − QC exclusion | 68 | 138 | blind record-level QC (§3) |
| **Analysed** | **259** | **516** | |

Every row was recomputed from `scripts/01_cohort.py` functions and `out/qc_designA.csv`, and it agrees with `out/gate1.csv`. The 5 cases lost between eligible and usable were all lost to the age criterion; none were lost to a missing file.

## 2. Included vs excluded: characteristics

Values are mean ± SD unless stated otherwise. Continuous variables were tested with Welch's t-test. Sex was tested with Pearson's chi-square without continuity correction.

| Group | Variable | Included | Excluded | p |
|---|---|---|---|---|
| All | n | 775 | 206 | |
| | Age, y | 68.3 ± 12.1 | 68.4 ± 12.8 | 0.914 |
| | Male, n (%) | 462 (59.6) | 129 (62.6) | 0.433 |
| | HR, bpm | 73.9 ± 12.3 | 73.9 ± 13.2 | 0.955 |
| | Beats, n | 12.2 ± 2.1 | 12.1 ± 2.0 | 0.695 |
| | P duration, ms | 145.7 ± 19.7 | 162.8 ± 28.1 | 1.2e−14 |
| | QRS, ms | 119.1 ± 26.6 | 162.2 ± 40.8 | 1.3e−34 |
| | QTcF, ms | 463.5 ± 55.2 | 498.7 ± 60.5 | 5.8e−13 |
| Cases | n | 259 | 68 | |
| | Age, y | 68.3 ± 11.9 | 68.2 ± 13.5 | 0.951 |
| | Male, n (%) | 152 (58.7) | 45 (66.2) | 0.261 |
| | HR, bpm | 77.2 ± 12.9 | 77.7 ± 14.7 | 0.796 |
| | Beats, n | 12.8 ± 2.2 | 12.6 ± 2.0 | 0.595 |
| | P duration, ms | 147.0 ± 19.2 | 160.4 ± 29.5 | 0.00064 |
| | QRS, ms | 128.5 ± 29.6 | 174.2 ± 42.0 | 7.1e−13 |
| | QTcF, ms | 485.4 ± 59.9 | 517.2 ± 64.5 | 0.00040 |
| Controls | n | 516 | 138 | |
| | Age, y | 68.2 ± 12.1 | 68.5 ± 12.6 | 0.857 |
| | Male, n (%) | 310 (60.1) | 84 (60.9) | 0.866 |
| | HR, bpm | 72.2 ± 11.7 | 72.1 ± 12.0 | 0.907 |
| | Beats, n | 11.9 ± 2.0 | 11.9 ± 2.0 | 0.942 |
| | P duration, ms | 145.0 ± 20.0 | 163.9 ± 27.5 | 2.1e−12 |
| | QRS, ms | 114.4 ± 23.6 | 156.3 ± 39.1 | 2.6e−24 |
| | QTcF, ms | 452.6 ± 49.2 | 489.5 ± 56.5 | 3.8e−11 |

Every one of the 981 records has P duration, QRS, QTcF, HR and beat count, because no record failed delineation or had a missing required feature. Three of these measurements — P duration, QRS and QT (behind QTcF) — are what QC checks. An excluded record either has one of them out of range or has QRS onset before P offset, which is a fiducial error that lengthens P and/or QRS. So the excluded-vs-included differences in P, QRS and QTcF partly come from how QC works. They should not be read as evidence that the excluded patients are physiologically different. Age, sex, HR and beat count are not QC inputs, and none of them differ.

## 3. Exclusions by QC reason (from the existing `qc_exclusions` logic)

*Flow* assigns each excluded record to its first failed rule, so the flow counts sum to the number excluded. *Any* counts every rule a record fails; a record can fail several rules, so these do not sum.

| Rule | All: flow | All: any | Cases: flow | Cases: any | Controls: flow | Controls: any |
|---|---|---|---|---|---|---|
| delineation failed | 0 | 0 | 0 | 0 | 0 | 0 |
| required feature / fiducial missing (6 rules) | 0 | 0 | 0 | 0 | 0 | 0 |
| qrs_onset_before_p_offset | 152 | 152 | 40 | 40 | 112 | 112 |
| range:hr | 0 | 0 | 0 | 0 | 0 | 0 |
| range:qrs_ms | 33 | 51 | 21 | 26 | 12 | 25 |
| range:qt_ms | 4 | 10 | 2 | 7 | 2 | 3 |
| range:p_ms | 17 | 26 | 5 | 10 | 12 | 16 |
| **Total excluded** | **206** | | **68** | | **138** | |

The script asserts that the flow counts sum to the excluded total in every group: 152 + 33 + 4 + 17 = 206; 40 + 21 + 2 + 5 = 68; 112 + 12 + 2 + 12 = 138. All counts match `out/qc_flow_designA.csv`.

## 4. Key check: does the exclusion rate differ between cases and controls?

| | Excluded / entered QC | Rate |
|---|---|---|
| Cases | 68 / 327 | 20.8 % |
| Controls | 138 / 654 | 21.1 % |

The rates do not differ: Pearson chi-square = 0.012, p = 0.912 (Fisher exact p = 0.934).

The **reasons** for exclusion do differ. Comparing the first-failed-rule distribution across the four non-empty rules gives chi-square = 17.7, df = 3, p = 0.0005. Two cells (range:qt_ms, 2 vs 2) have small expected counts, so treat this p as approximate. Cases are excluded relatively more often for QRS outside 50–200 ms (21/68 vs 12/138 as the first failure). Controls are excluded relatively more often for QRS onset before P offset (112/138 vs 40/68). `qc_exclusions.md` already reports this.

## Summary (post hoc)

1. QC removed the same proportion of cases and controls (20.8 % vs 21.1 %, p = 0.91). Within each group, excluded and included records do not differ in age, sex, HR or beat count.
2. Excluded records have longer P, QRS and QTcF. This is largely by construction, because those measurements, or the fiducials behind them, are what QC tests.
3. The mix of exclusion reasons differs by group (p = 0.0005). Wide-QRS exclusions fall more on cases and P/QRS-overlap exclusions more on controls, so the QRS comparison after QC may be conservative for cases.
