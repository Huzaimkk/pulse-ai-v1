# PULSE-AI v1

Automated, validated extraction of named atrial and ventricular markers from the
10-second 12-lead ECG (P duration, PTFV1, P axis, interatrial block, QRS, QT/QTc,
Tpeak–Tend, ST level, heart rate, RMSSD).

Code for the paper *"PULSE-AI v1: Validated Extraction of Atrial and Ventricular
Markers from the 10-Second 12-Lead ECG, Tested on Left Atrial Abnormality in PTB-XL"*
(Department of Electronic Engineering, Dawood University of Engineering and
Technology, Karachi).

## Setup

```
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt     # Windows
```

Python 3.11. Data are not included: download PTB-XL 1.0.3, LUDB 1.0.1 and QTDB
1.0.0 from PhysioNet into `data/` (paths in `config.yaml`).

## Layout

| Folder | Contents |
|---|---|
| `pulse_ecg/` | filtering, median beat, wavelet delineation, marker registry, features, quality control |
| `scripts/` | cohort building and delineation validation |
| `paper/notes/` | analysis plan, analysis scripts and their written results |
| `out/` | result tables (CSV) produced by the scripts |
| `figures/` | figure scripts and figures |
| `paper/` | LaTeX source; `make_numbers.py` writes every number in the text from `out/` |
| `tests/` | unit tests on synthetic ECGs with known answers |

All tunable constants live in `config.yaml`. Run the tests with `pytest`.

## Data

PTB-XL, LUDB and QTDB are open-access databases on PhysioNet. Please cite them
as given in the paper.
