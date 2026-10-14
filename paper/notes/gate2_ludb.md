# Gate 2 — LUDB delineation validation

> **Decision: `p_span_ms` = 100 (column A) is the chosen value and is what
> `config.yaml` and `out/gate2_*.csv` now hold.** 150 (column B) was tested and
> rejected: on held-out LUDB sinus records it moved the P-duration bias from
> +5.7 to +34.0 ms and κ at 120 ms from 0.26 to 0.01. An accurate measurement
> with documented missing data is usable; a biased measurement on every record
> is not. Column B is kept below only as the record of that test — **it is not
> current**. The cost of choosing 100 is measured in `p_missingness.md`.

Held-out LUDB half (split seed 42), **sinus records only**: 87 records, 824
reference beats, 711–712 scored P waves. Produced by
`scripts/02_validate_delineation.py`; current outputs in `out/gate2_*.csv`.
Delineator code identical in both columns; only `wavedet.p_span_ms` differs.

**Column A was measured at `p_span_ms` = 100 (current config.yaml). Column B
was measured at `p_span_ms` = 150 (tested, rejected).** The column-A numbers were
re-run at 100 with the current code and reproduce the earlier `out/` files
exactly, so every difference below is due to `p_span_ms` alone. `out/gate2_*.csv`
were regenerated at 100 on 2026-09-29 and match column A.

## Primary reference: k-of-N rule (same rule the pipeline measures with)

| Metric | A: p_span 100 | B: p_span 150 | CSE 2-SD tol. | Martínez 2004 (QTDB) |
|---|---:|---:|---:|---:|
| Within-record SD, P onset (ms) | 15.2 | **31.2** | 10.2 | 14.8 |
| Within-record SD, P offset (ms) | 12.2 | 13.5 | 12.7 | 12.8 |
| Within-record SD, QRS onset (ms) | 12.0 | 12.0 | 6.5 | 7.7 |
| Within-record SD, QRS offset (ms) | 6.4 | 6.4 | 11.6 | 8.7 |
| Within-record SD, T offset (ms) | 24.5 | 24.5 | 30.6 | 18.1 |
| Mean P-onset error (ms) | −4.6 | **−29.9** | | |
| P-duration bias, per-record median, mean ± SD (ms) | +5.7 ± 23.4 | **+34.0 ± 35.1** | | |
| P-duration MAE, per record (ms) | 16.3 | **38.1** | | |
| Records measuring P ≥ 120 ms (delineator) | 70 % (61/87) | **87 % (76/87)** | | |
| Records with P ≥ 120 ms (reference) | 66 % | 66 % | | |
| κ at 120 ms | 0.26 | **0.01** | | |

## Secondary reference: earliest onset / latest offset over annotated leads

| Metric | A: p_span 100 | B: p_span 150 |
|---|---:|---:|
| Within-record SD, P onset (ms) | 19.8 | **33.0** |
| Within-record SD, P offset (ms) | 13.1 | 13.9 |
| Within-record SD, QRS onset (ms) | 13.2 | 13.2 |
| Within-record SD, QRS offset (ms) | 7.8 | 7.8 |
| Within-record SD, T offset (ms) | 30.8 | 30.8 |
| P-duration bias, per-record median, mean ± SD (ms) | −8.0 ± 28.8 | **+20.4 ± 37.9** |
| Records measuring P ≥ 120 ms (delineator) | 70 % | **87 %** |
| Records with P ≥ 120 ms (reference) | 80 % | 80 % |
| κ at 120 ms | 0.24 | 0.16 |

## Status

Widening `p_span_ms` from 100 to 150 (needed so synthetic P waves of
155–200 ms are delineated at all) **makes P onset and P duration substantially
worse on real ECGs**. P onset moves on average ~25 ms earlier: the wider span
admits a slope that precedes the true P onset. QRS and T boundaries are
unchanged. Per instruction, `p_span_ms` was not tuned further. The trade-off
between wide-P detection and P-onset accuracy on LUDB is unresolved and needs a
decision before these numbers go into the paper. The Gate 2 numbers in
`methods_delineation_and_endpoints.md` are the column-A (p_span 100) values.

**Resolved (2026-09-29):** 100 chosen; see the decision box at the top and
`p_missingness.md` for its cost on the Design A cohort.
