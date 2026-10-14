# POST-HOC revision: methods facts, read from code

**Revision material, not pre-specified.** Every fact below comes from reading the code at the cited `file:line`. Points 5 and 7 also add measurements, and those are labelled as such. No code was changed.

## 1. Synthetic ECGs (`pulse_ecg/synth.py`)

**Model.** This is a deterministic, non-physiological template, not a dynamical model such as ECGSYN.

- **Wave shapes.** Every wave is a half-sine, or a chain of half-sines, with exact support, so onsets and offsets are known to the sample (`synth.py:4-8`, `synth.py:20-24`).
- **Limb leads.** Leads I and II are projections of a frontal-plane axis at 0° and 60°. III, aVR, aVL and aVF are derived from them by Einthoven's and Goldberger's relations, so III = II − I holds exactly and the P axis has a known true value (`synth.py:27-41`).
- **Precordial leads.** These use fixed gains per wave (P, QRS, T). V1 has a negative QRS, giving an rS complex (`synth.py:15-17`).
- **QRS.** Built from Q, R and S half-sines covering 15 %, 50 % and 35 % of the QRS width, with amplitudes −0.10, 1.0 and −0.25 × R (`synth.py:91-94`).
- **T wave.** One half-sine of width 0.5 × (QT − QRS), ending at QRS onset + QT (`synth.py:68`, `synth.py:95`).
- **PTFV1 truth.** With `biphasic_v1=True`, V1's P wave is positive and then negative. The known answer is PTFV1 = −(terminal amplitude in µV) × (terminal duration in ms) (`synth.py:85-88`, `synth.py:115-116`).
- **Noise.** Optional white Gaussian noise, `noise_mv` (default 0), seeded (`synth.py:97-98`).
- **Rhythm.** Strictly regular: beat k sits at a fixed RR = 60/HR, so there is no RR variability and no ectopy (`synth.py:66`, `synth.py:70-75`).

**Parameters and defaults** (`synth.py:44-49`):

| Parameter | Default |
|---|---|
| hr | 60 bpm |
| duration_s | 10 s |
| fs | 500 Hz |
| p_duration_ms | 100 |
| pr_ms | 160 |
| qrs_duration_ms | 90 |
| qt_ms | 400 |
| p_amp_mv | 0.15 |
| r_amp_mv | 1.0 |
| t_amp_mv | 0.3 |
| p_axis_deg | 50 |
| qrs_axis_deg | 45 |
| t_axis_deg | 45 |
| v1_terminal_frac | 0.5 |
| v1_terminal_amp_mv | 0.1 |
| noise_mv | 0 |

The function returns the true fiducials and feature values with `return_truth=True` (`synth.py:106-117`).

**What they were used for.** Only for tests with known answers and for characterising the delineator. They were never used in training or in any cohort analysis.

- **Unit and integration tests.** There are 35 `synth(` calls across 8 files under `tests/`: `test_atrial.py` (13), `test_synth.py` (7), `test_ventricular.py` (4), `test_extract.py` (3), `test_rate.py` (3), `test_delineate.py` (2), `test_io_preprocess.py` (2) and `test_features.py` (1). One noisy variant uses `noise_mv=0.02` (`tests/test_delineate.py:18`).
- **Delineator characterisation notes:**
  - P-duration compression, "measured = 0.782 × true + 27.9 ms" (`paper/notes/p_duration_compression.md:3-6`, cited in `gate3_sanity.md:68`)
  - P-axis flip (`paper/notes/p_axis_flip.md:1`)
  - P missingness for P ≥ 155 ms at `p_span_ms` 100 (`paper/notes/p_missingness.md:27`)
  - the `p_span_ms` choice (`paper/notes/gate2_ludb.md:54`)

## 2. Ectopic / premature beat handling before the median beat

**There is no RR-based rejection of premature or ectopic beats, and no prematurity threshold.**

- **R-peaks** (`delineate.py:70-100`). These come from the integrated energy of the multilead derivative. The detector applies only a refractory distance of 250 ms and an amplitude gate of 0.3 × the 90th-percentile peak height (`config.yaml:31`; `delineate.py:78-83`). It does not classify beats.
- **Median beat** (`delineate.py:103-115`). Every detected R-peak whose −0.40 s / +0.50 s window fits in the record is used (`config.yaml:32`; `delineate.py:107`). The template is the sample-wise median of all those beats (`delineate.py:111`). At least `min_beats: 3` complete beats are required (`config.yaml:32`).
- **Template correlation.** A beat-to-template correlation is computed (`delineate.py:112-114`) and stored as `beat_template_corr` (`delineate.py:44`, `delineate.py:65`). Nothing reads it: it is not used to reject beats or records, and no other module in `pulse_ecg/` references it.
- **Wavedet fiducials.** Per-beat marks are combined per record with a plain `nanmedian` of (mark − R) over all beats (`delineate.py:143-145`). The preceding RR only sizes each beat's search window. For P, the window starts at R − min(0.45 × preceding RR, 350 ms) (`wavedet.py:125-130`, `wavedet.py:175`; `config.yaml:59-60`). It does not exclude any beat.
- **RR-based features.** `mean_hr_bpm`, `rmssd_ms` and the QTc RR use every RR interval, unfiltered (`features/rate.py:25-26`, `rate.py:112-116`, `rate.py:140-141`; `features/ventricular.py:46-50`).
- **Practical effect.** The median gives robustness only while ectopic beats are a minority of a record's beats. The analysed records are all PTB-XL "SR" records, but SR does not exclude occasional ectopic beats.

## 3. QTDB handling (Gate 2, secondary)

- **Sampling rate.** QTDB is not resampled; it is processed at its native 250 Hz.
  - Gate 2 tags these rows `fs=250` (`scripts/02_validate_delineation.py:168`).
  - The band-pass and notch are designed for each record's own fs (`preprocess.py:13-24`, `preprocess.py:30-31`).
  - The wavelet transform is defined in the 250 Hz numbering of Martínez 2004 and accepts fs = 250 × 2ⁿ, shifting scales for 500 Hz (`wavedet.py:49-54`). So 250 Hz is wavedet's native rate, and nothing in the method is adapted beyond fs-based ms → sample conversions.
  - The errors are therefore quantised to 4 ms, against 2 ms for LUDB (`02_validate_delineation.py:14-15`).
- **Leads.** Both of the record's two leads, whatever their names, enter the multilead rule (`io.py:94-95`). The annotation is not per lead (`io.py:89-90`).
- **Annotations.** `q1c`, cardiologist 1's manual annotation (`io.py:88-96`; `02_validate_delineation.py:11`). The annotated span (~30 beats) is cut into 10 s windows (`02_validate_delineation.py:156-166`). A beat is scored only if its whole P-QRS-T sits inside the window (`02_validate_delineation.py:165-166`).
- **Records.** `RECORDS` lists 105, and all 105 have a `.q1c` file. 104 were scored (3,178 beats, per `out/gate2_beats.csv`). `sel232` was skipped because its q1c has no annotated QRS (`02_validate_delineation.py:153-155`). QTDB is labelled `half="external"`, and no constants were tuned on it (`02_validate_delineation.py:168`).

## 4. Design B: "AF elsewhere"

**Definition**

- **Cases** are sinus records with no AFIB code of their own, from patients who have **AFIB in at least one other PTB-XL record**. AFIB is the SCP code `AFIB` (`config.yaml:163`) at any likelihood under `min_likelihood: 0` (`scripts/01_cohort.py:49`, `01_cohort.py:67-71`).
- **Controls** are sinus records of patients who never carry AFIB in any record (`01_cohort.py:76`).
- Selection then follows Design A: one record per patient, the earliest eligible (`01_cohort.py:76`); age ≥ 40; 1:2 matching.

**Timing.** The timing relative to the sinus record is knowable but is not used. The code supports `before_first_af`, which keeps only sinus records dated before the patient's first AFIB record (`01_cohort.py:72-73`). The configured setting is `af_record_timing: any` (`config.yaml:164`), so the sinus record may come before or after the AF record.

**Measured post hoc** from `out/cohort.csv` and `recording_date`, for the 118 matched Design B cases:

| Sinus record relative to first AFIB record | Cases |
|---|---|
| Before | 51 |
| After | 67 |
| Same date | 0 |

So 57 % of Design B "cases" were recorded after AF had already been documented. For those records Design B is not a prediction design.

## 5. Notch filter

**Settings.** The notch is set at 50 Hz with Q = 30 (−3 dB bandwidth 50/30 ≈ 1.7 Hz), applied together with a 0.5–40 Hz Butterworth band-pass of order 4 (`config.yaml:30`).

- The notch is `scipy.signal.iirnotch`, applied zero-phase with `sosfiltfilt` before the band-pass (`preprocess.py:21-23`, `preprocess.py:33-35`).
- The notch is skipped only if 50 Hz ≥ Nyquist (`preprocess.py:21`).
- The same filter is used for every source: PTB-XL, LUDB, QTDB and the device (`preprocess.py:1-2`).

**Does PTB-XL need it?** PTB-XL was recorded in Germany, where mains is 50 Hz, so 50 Hz is the correct frequency. Two post-hoc measurements show the notch has almost nothing to remove:

- **The band-pass already attenuates 50 Hz.** Its response at 50 Hz is 0.358 per pass, which is −17.9 dB after zero-phase filtering (both passes). At 60 Hz it is −30.8 dB.
- **The raw records carry almost no 50 Hz content.** Across the 775 analysed PTB-XL records (Welch PSD, lead-averaged), the 50 Hz peak sits a median 0.5 dB (IQR 0.1–1.0) above the 45–55 Hz floor. It exceeds 10 dB in 6/775 (0.8 %) records and 20 dB in none.

Measured by a scratch script in the session scratchpad, not part of the repo.

## 6. QTc abnormal threshold

**There is a single threshold, ≥ 460 ms, for both sexes.** It is not sex-specific.

- `config.yaml:102`: `qtc_ms: {normal: [390, 450], abnormal: [">=", 460]}`, cited to rautaharju2009aha4.
- `config.yaml:100-101` records the choice: "QTc prolonged >= 450 ms (men) / >= 460 ms (women) … One threshold per feature: 460 ms is prolonged for both sexes."
- The same entry is attached to `qtc_bazett_ms` (`features/ventricular.py:102-108`) and `qtc_fridericia_ms` (`ventricular.py:111-117`) via `QTC = V["qtc_ms"]` (`ventricular.py:26`).
- Both features use mean RR over all beats (`ventricular.py:46-50`).
- QTc enters the analysis as a continuous value. No QTc threshold call appears in `out/group_threshold_calls_designA.csv`.

## 7. Python and package versions

Python 3.11.16 in `.venv`. The venv has no pip, so `paper/versions.txt` was written with `uv pip freeze --python .venv/Scripts/python.exe` (49 packages). Its output matched an `importlib.metadata` listing exactly.

Key versions:

| Package | Version |
|---|---|
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| pandas | 2.3.3 |
| scikit-learn | 1.9.1 |
| wfdb | 4.3.1 |
| joblib | 1.6.0 |
| matplotlib | 3.11.2 |
| PyYAML | 6.0.3 |
