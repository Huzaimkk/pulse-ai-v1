# P-axis 180° flip on synthetic 110 and 130 ms P waves

Neither the fixture nor `atrial.py`: the multilead wavedet puts the global QRS onset
inside the P wave (qrs_onset 130/126 samples < p_offset 154/153; true ≈ 177). Then
`delineate.py` re-references the median beat to "the isoelectric run before QRS
onset", so the baseline lands on the P wave (+0.05 mV in lead I). The P areas in I
and aVF both go negative, and `atan2` correctly returns 50° − 180° = −130°. It was
not fixed, per instruction (delineator untouched). Found 2026-09-29, `p_span_ms` 100.
