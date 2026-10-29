# Multistate lifestyle analysis

This repository contains R scripts for analysing lifestyle, social support and trajectories of prostate disease and depression in UK Biobank. The scripts are organised by analysis module and document the workflow and key parameters used in the study.

The multistate model includes Healthy, benign prostatic hyperplasia (BPH), Depression, prostate cancer (PCa) and Death.

## Analysis scripts

All scripts are located in `scripts/`.

| File | Analysis module |
|---|---|
| `00_gp_screening.R` | GP code filtering, disease-date matching and earliest-date aggregation by participant. |
| `01_baseline_first_occurrence_GP.R` | Baseline cohort preparation and integration of first-occurrence, death and GP records. |
| `02_inpatient_dates.R` | Disease-date extraction from inpatient records. |
| `03_cancer_registry_dates.R` | Prostate cancer date extraction from cancer registry records. |
| `04_self_report_dates.R` | Disease-date extraction from self-reported diagnoses. |
| `05_followup_and_exclusions.R` | Outcome-date integration, follow-up construction and cohort exclusions. |
| `06_covariates_imputation.R` | Baseline covariate preparation, chained-equation imputation and risk-score construction. |
| `07_paths_and_risk_sets.R` | Multistate pathway and transition-specific risk-set construction. |
| `08a_individual_factor_cox.R` | Transition-specific Cox models for individual lifestyle-related factors. |
| `08b_categorical_score_cox.R` | Cox models for categories of the cumulative lifestyle risk score. |
| `08c_continuous_score_cox.R` | Cox models per additional higher-risk factor. |
| `08d_component_sensitivity_cox.R` | Dietary and social-support component sensitivity analyses. |
| `09_descriptive_tables.R` | Baseline characteristics, incidence rates, state durations and annual summaries. |
| `10_state_occupancy.R` | State occupancy probabilities by age and lifestyle risk-score strata. |
| `11_scenarios_and_bootstrap.R` | Hypothetical exposure scenarios, predicted years gained and participant-level bootstrap. |
| `11b_individual_predictions.R` | Individual predicted gains and an illustrative high-benefit profile. |
| `12_timing_sensitivity.R` | Sensitivity analyses of same-day events, landmarks and diagnostic intervals. |
| `13_smoking_implementation.R` | Partial-implementation sensitivity analysis for the smoking-risk indicator. |
| `14_plot_figures.R` | Optional plotting from saved statistical results. |

## Data and usage

UK Biobank participant data, derived datasets, fitted models and intermediate objects are not included in this repository. By default, scripts read inputs from `data/raw/` and save outputs to `results/`. Local paths, input structures and software settings may require adaptation before analysis.

Set the working directory to the repository root. Run the scripts individually in their numbered order, completing the required upstream steps first. Sections can also be run step by step in RStudio to inspect intermediate objects.

## Software

Recorded versions: R 4.4.1; `data.table` 1.16.2; `survival` 3.6-4; `mstate` 0.3.3; `mice` 3.16.0; and `ggplot2` 3.5.1 for optional plotting.

These versions are provided for reference, not as a newly validated environment. A previously observed ABI compatibility warning between `lme4` 1.1-35.5 and `Matrix` 1.6-5 remains unresolved. The final code packaging did not rerun the study analyses.
