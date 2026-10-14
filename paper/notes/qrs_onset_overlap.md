# QRS onset placed before P offset on the median beat (Design A)

**Conclusion: not a synthetic-fixture artifact. It affects 15.5 % of Design A records
(152/981), and QRS duration, PR interval and P axis are all compromised in them.**

Measured 2026-09-29 with the current `config.yaml` (delineation.method = wavedet,
PTB-XL records500, fs 500 Hz), production `context()` + `REGISTRY.compute`, script
`paper/notes/qrs_onset_overlap.py`, per-record output
`out/qrs_onset_overlap_designA.csv`. No errors; p_offset and qrs_onset present on
all 981 records. Nothing in the delineator, `atrial.py`, `ventricular.py` or config
was changed.

## Frequency of qrs_onset < p_offset

| group   | n   | overlap | %      |
|---------|-----|---------|--------|
| all     | 981 | 152     | 15.5 % |
| case    | 327 | 40      | 12.2 % |
| control | 654 | 112     | 17.1 % |

A further 20 records have qrs_onset == p_offset exactly (gap 0).

## Gap qrs_onset − p_offset (ms)

| group   | median | p90  | p10  |
|---------|--------|------|------|
| all     | 37     | 69   | −6   |
| case    | 42     | 75   | −4   |
| control | 34     | 66   | −8   |

Range −70 to +144 ms. The boundaries normally sit close: 103/981 records have a gap
of 0–10 ms. Depth of the overlap in the 152 affected records: 0–5 ms 51, 5–10 ms 25,
10–20 ms 31, 20–30 ms 17, 30–70 ms 27.

## Effect on the affected records (medians, overlap vs no overlap)

| feature                         | overlap (n=152) | no overlap (n=829) |
|---------------------------------|-----------------|--------------------|
| qrs_duration_ms                 | 157             | 114                |
| PR, qrs_onset − p_onset (ms)    | 140             | 188                |
| p_duration_ms                   | 158             | 144                |
| records with \|p_axis_deg\| > 120° | 60 (39 %)       | 16 (2 %)           |

QRS onset is pulled early, so QRS duration is inflated and PR is shortened.
The P axis flips on real records: 60 of the 152 affected records have
|p_axis_deg| > 120° (all 60 are < −90° or > 150°), against 16 of 829 unaffected
records. The flip tracks the depth of the overlap: 50 of the 68 records with an
overlap ≥ 15 ms flip, against 10 of the 84 with a shallower overlap. This is the
same mechanism as the synthetic 110/130 ms case in `p_axis_flip.md`: the
re-reference window before qrs_onset falls on the P wave. Example flipped records
(ecg_id, gap ms, p_axis_deg): 21206 −27 −158; 8050 −40 −125; 14292 −44 −135;
3265 −37 −125; 541 −38 −170; 19440 −52 +180.

## Case/control imbalance

The overlap is more frequent in controls (17.1 %) than cases (12.2 %), and so are
flipped axes (48 controls vs 12 cases among affected records). Any case–control
comparison of QRS duration, PR or P axis on Design A is therefore biased by
delineation error, not just noisier.

The per-record values for all 152 affected records (p_offset, qrs_onset, gap, PR,
p_axis_deg, p_duration_ms, qrs_duration_ms) are in the CSV above (filter gap_ms < 0).
