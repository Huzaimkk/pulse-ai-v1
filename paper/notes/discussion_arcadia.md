# Discussion — ARCADIA paragraph (draft)

Goes in Section IX (Discussion), after the clinical-meaning paragraph.
Numbers verified against the PubMed abstract (PMID 38324415) on 2026-09-28.

## Facts from the trial (Kamel et al., JAMA 2024;331:573–581, `kamel2024arcadia`)

- Population: 1015 patients with **cryptogenic stroke** and atrial cardiopathy,
  no AF at randomisation. It was a **secondary**-prevention trial.
- Atrial cardiopathy was defined by any one of these: **PTFV1 > 5000 µV·ms in V1**,
  NT-proBNP > 250 pg/mL, or LA diameter index ≥ 3 cm/m².
- Apixaban vs aspirin. Recurrent stroke was 40 vs 40 (4.4 %/yr in each arm),
  **HR 1.00 (95 % CI 0.64–1.55)**. The trial was **stopped for futility** at a
  planned interim analysis.

## The argument

ARCADIA asked a treatment question: *given* atrial cardiopathy, does
anticoagulation beat aspirin? It did not ask whether the markers identify
elevated risk. Every participant had already had a stroke, and both arms had
the same high recurrence rate. The marker–risk association this paper relies on
comes from cohort evidence (BAYES, `martinezselles2020bayes`), and ARCADIA does
not test it. So the trial does not weaken the case for measuring these markers.
What it removes is any claim about **what to do** once they are found.

## Draft text

> The ARCADIA trial [kamel2024arcadia] randomised 1015 patients with cryptogenic
> stroke and atrial cardiopathy, defined in part by PTFV1 > 5000 µV·ms, to
> apixaban or aspirin, and stopped for futility with identical recurrence rates
> (HR 1.00, 95 % CI 0.64–1.55). ARCADIA tested a treatment, not a marker: it
> shows that anticoagulation did not improve outcomes in this population, not
> that atrial cardiopathy fails to identify risk, which rests on independent
> cohort evidence [martinezselles2020bayes]. The present work makes no treatment
> claim. Its contribution is measurement: whether these markers can be extracted
> automatically and reliably, which is a precondition for any future trial that
> selects patients by them.

## Rule (also in CLAUDE_1.md)

Never write that patients with these markers should be anticoagulated, or imply
it. That includes softer versions such as "may benefit from anticoagulation" and
"candidates for anticoagulation".

## Open item for Fig. 3

Stage 3's body reads *"preventable with anticoagulation / only in patients
identified"*. It refers to AF, where anticoagulation is guideline-indicated.
But sitting in a figure about Stage-1 markers, it can be read as recommending
anticoagulation for patients identified by those markers. Consider rewording,
for example to *"anticoagulation indicated once AF is documented"*.
