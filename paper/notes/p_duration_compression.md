# P-duration compression of the delineator

Measured on noise-free synthetic 12-lead ECGs (`pulse_ecg.synth`, half-sine P,
HR 60, PR = max(160, P + 50) ms), 500 Hz, through the full pipeline
(`extract` → multilead wavedet, `p_span_ms: 100`, the chosen value; re-measured
2026-09-29). Measured value is the reported `p_duration_ms` endpoint.

| True P (ms) | Measured P (ms) | Error (ms) |
|---:|---:|---:|
| 80  | 92  | +12 |
| 90  | 98  | +8  |
| 100 | 106 | +6  |
| 110 | 112 | +2  |
| 130 | 130 | 0   |
| 150 | 146 | −4  |

170 and 200 ms are not listed: at `p_span_ms` 100 a clean half-sine P ≥ 155 ms is
not delineated (NaN); on real records this never happened (0/981, `p_missingness.md`).

Least-squares fit (80–150 ms): measured = 0.782 × true + 27.9 ms (r = 0.998).
Slope **0.782**; the identity line is crossed at about 128 ms.

The compression attenuates the case-control difference, so the reported effect size is conservative.

The fit has an intercept of **+27.9 ms** as well as a slope of 0.782, and the two
consequences point in opposite directions:

- **Continuous case-control comparison (primary endpoint).** Only the slope
  matters for a difference between groups: a true difference Δ is measured as
  0.782 Δ. The compression attenuates the difference, so the reported effect
  size is conservative.
- **The ≥ 120 ms binary threshold (IAB).** Here the intercept matters, and it
  works the other way. Below the ~128 ms crossover the delineator reads long: a
  true 118 ms P wave measures 0.782 × 118 + 27.9 = 120 ms. The threshold
  therefore fires early, and IAB prevalence is overestimated.

Cross-reference to LUDB (Gate 2, `methods_delineation_and_endpoints.md`, and
`gate2_ludb.md`): on held-out sinus records the k-of-N *reference* built from the
cardiologists' annotations already puts 66 % of records at ≥ 120 ms, against a
plausible prevalence of 10–20 %. That figure involves no delineator output, so it
belongs to the annotations and the reference rule. The delineator's *measured*
prevalence on the same records is higher still (70 % at `p_span_ms` 100; see
`gate2_ludb.md` for 150). Its positive bias at short widths adds to the
reference's inflation rather than explaining it. The delineator's bias and the
reference's inflation both push the 120 ms dichotomy upward. This is a further
reason `p_duration_ms` (continuous) is the primary endpoint and `iab_class` is
descriptive only.
