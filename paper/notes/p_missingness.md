# P-wave feature missingness at `p_span_ms` = 100 (Design A)

**For Limitations.** Measured 2026-09-29 with `paper/notes/p_missingness.py`: the
production `extract()` (preprocess → multilead wavedet → feature registry), current
`config.yaml` (`wavedet.p_span_ms: 100`), on every Design A record in
`out/cohort.csv`. Per-record output, including the NaN reason for each feature, is in
`out/p_missingness_designA.csv`. No values were imputed and no records were dropped
before counting.

## Result

| Feature | Overall NaN | Cases NaN | Controls NaN |
|---|---:|---:|---:|
| `p_duration_ms` (primary endpoint) | 0 / 981 (0.0 %) | 0 / 327 (0.0 %) | 0 / 654 (0.0 %) |
| `ptfv1_uv_ms` | 1 / 981 (0.1 %) | 1 / 327 (0.3 %) | 0 / 654 (0.0 %) |
| `p_axis_deg` | 0 / 981 (0.0 %) | 0 / 327 (0.0 %) | 0 / 654 (0.0 %) |

The single missing PTFV1 is case ecg_id 15703 (reason: "no single-lead P boundaries
in V1"). Its global P duration (112 ms) and P axis were measured.

**No differential missingness.** The case rate is not higher than the control rate
for any of the three features in any way that matters (one record). Missing data on
the primary endpoint is not a limitation at `p_span_ms` = 100 in this cohort.

## Why this differs from what the synthetic test predicted

With `p_span_ms` = 100, noise-free synthetic half-sine P waves of ≥ 155 ms are *not
delineated* (checked again today: 150 ms → 146 ms measured; 155, 170 and 200 ms →
NaN, "requires p_onset, p_offset, not delineated"). That was the reason 150 was
tried. On real PTB-XL records the failure does not show up as NaN: every Design A
record gets a global P onset and offset. Measured values go up to 234 ms
(cases: median 147, 5th–95th pct 117–194 ms; controls: median 144, 118–192 ms).
Real P waves are not clean two-slope half-sines, so the delineator always finds a
qualifying slope pair within 100 ms.

## What the Limitations section should therefore say

The cost of choosing 100 is **not** missing data. The possible cost is a
*misplaced* boundary. For a very wide real P wave whose outer slopes are more
than ~100 ms apart, the delineator may place the onset or offset on an inner slope,
which gives a shorter value instead of NaN. This cohort cannot measure that because
it has no reference annotation. The evidence on it comes from held-out LUDB (Gate 2,
`gate2_ludb.md`, column A): the P-duration bias is +5.7 ± 23.4 ms, so there is no net
shortening, but LUDB contains few P waves ≥ 155 ms, so that range is weakly tested.
Suggested wording: *"P-wave boundaries were located with a 100 ms slope-span
constraint, chosen because a wider span biased P duration by +34 ms on LUDB. No
Design A record lacked a P-duration measurement (0/981; PTFV1 1/981), so missing
data did not differ between cases and controls. On synthetic signals this constraint
fails to delineate P waves ≥ 155 ms, and on real records it may truncate rather than
omit such waves. The LUDB validation contains few P waves in this range."*

`p_duration_compression.md` was re-measured at `p_span_ms` = 100 on 2026-09-29.
