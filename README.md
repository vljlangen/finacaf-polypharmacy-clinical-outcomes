# FinACAF Polypharmacy Clinical Outcomes

Nationwide retrospective cohort analysis of incident atrial fibrillation (FinACAF, 2007-2018) evaluating associations between medication count and risks of ischemic stroke, mortality, intracranial hemorrhage, and major bleeding.

## Repository contents
- `Polypharmacy_outcomes_260426.R`: primary analysis and figure generation.
- `km_panel_polypharmacy_clinical_outcomes.R`: Kaplan-Meier panel figure script.

## Data
Original register data are not included. Scripts require external source files provided by data permit holders.

## Run
Set input file environment variables, then run the R scripts:
- `POLYPHARMACY_MAIN_COHORT_SAV`
- `POLYPHARMACY_MEDICATION_SAV`
- `POLYPHARMACY_POST_AF_CATEGORY_SAV`
