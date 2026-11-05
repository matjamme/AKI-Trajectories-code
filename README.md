# AKI Trajectories — analysis code

Code for *Competing-risk Trajectory Phenotyping of Acute Kidney Injury in Critical Illness*.

## Scripts

| Script | Analysis |
| --- | --- |
| `01_JLCM_fitting.R` | Competing-risk quadratic Jointlcmm model, SCr winsorisation at the 99th percentile, model selection |
| `02_External_class_assignment.R` | Internal evaluation and external posterior class assignment in MIMIC-IV and Genevadb |
| `03_XGBoost_training.R` | XGBoost model development (M0–M4 and SCr-only) |
| `04_XGBoost_external_validation.R` | Landmark-aligned training and external performance evaluation in MIMIC-IV |
| `05_Landmark_Cox_bootstrap.R` | Cox regression and paired bootstrap C-index analysis |
| `06_Enrichment_analysis.R` | Prognostic enrichment and sample-size calculations |

## Input data

The source datasets are accessed independently under their respective access terms. No patient-level data are distributed with the code.

- [MIMIC-IV source and derived-concept code](https://github.com/MIT-LCP/mimic-code)
- [AmsterdamUMCdb resources](https://github.com/AmsterdamUMC/AmsterdamUMCdb)

Input datasets must contain the study-specific AKI-SCr cohorts and the variables expected by each script.

## JLCM fitting

```bash
Rscript 01_JLCM_fitting.R --db_path=/secure/AUMC_datA_AKI_SCr.rds --dataset=AUMC --out_dir=/secure/models --g_min=2 --g_max=5 --horizon_d=14 --time_var=time_s
```

## External class assignment

The evaluation entry point dispatches to cohort-specific implementations in `internal/`:

```bash
Rscript 02_External_class_assignment.R --cohort=AUMC --aumc=/secure/AUMC.rds --models=/secure/models --pattern=G3
Rscript 02_External_class_assignment.R --cohort=MIMIC --aumc=/secure/AUMC.rds --external=/secure/MIMIC.rds --models=/secure/models --out=/secure/results/mimic --pattern=G3
Rscript 02_External_class_assignment.R --cohort=GENEVA --aumc=/secure/AUMC.rds --external=/secure/Geneva.rds --models=/secure/models --out=/secure/results/geneva --pattern=G3
```

## XGBoost

Set `AKI_BASE_DIR` to the secured project folder containing `DF/` and `Modele Joint/`:

```bash
AKI_BASE_DIR=/secure/TrajAKI_Ventil Rscript 03_XGBoost_training.R
AKI_BASE_DIR=/secure/TrajAKI_Ventil Rscript 04_XGBoost_external_validation.R
```

The cumulative observation windows run from 24 to 192 hours after ICU admission. In the evaluation script, patients must be alive and still in the ICU at the end of each window. The original internal day codes `0` through `7` correspond to 24 through 192 hours.

The XGBoost pipeline uses five-fold cross-validation, early stopping, and bootstrap intervals for discrimination and calibration metrics.

## Landmark Cox models and paired bootstrap

`05_Landmark_Cox_bootstrap.R` reads the frozen M2 predictions, clinical covariates, MIMIC-IV ICU admission/discharge/death timestamps, and SCr time series. At 24, 48 and 72 hours, it estimates the C-index of the clinical model alone and after adding the two independent class probabilities `pC1` and `pC2`. At 48 hours, it also adds admission SCr and the change through 48 hours to the same clinical model, using the same analysis population. Models are re-fitted in each of 1,000 paired patient-level bootstrap samples. Separately, the retrospective clinical reference model is fitted in AmsterdamUMCdb and externally evaluated in MIMIC-IV using ICU mortality; its external C-index interval resamples MIMIC-IV stays with the AmsterdamUMCdb model held fixed.

```bash
Rscript 05_Landmark_Cox_bootstrap.R \
  --pred=/secure/LANDMARK_pred_dt.csv \
  --clinical=/secure/comparison_variables_MIMIC.rds \
  --baseline=/secure/mimic31_aki_admission_base_v3.rds \
  --scr=/secure/timeseries_mimic.rds \
  --aumc_clinical=/secure/comparison_variables_AUMC.rds \
  --aumc_baseline=/secure/aumc_bf.rds \
  --out=/secure/results/cox \
  --boot=1000
```

Outputs include `Table_S19_landmark_clinical_vs_M2.csv`, `Table_S19_retrospective_clinical_reference.csv`, separate bootstrap replicate files for each landmark, and the session information.

## Mortality enrichment and trial sample-size calculations

`06_Enrichment_analysis.R` calculates subsequent ICU mortality enrichment at the 24- and 48-hour landmarks. Its first module evaluates M2 top-10% and top-20% selection, enrichment factors, and projected two-arm trial sample sizes over a range of relative risk reductions. Its second module compares M2 selection with RRT initiated by the landmark and maximum SCr-based KDIGO stage 3, including equal-selection-fraction comparisons and analyses among patients without RRT. Both modules use 2,000 patient-level bootstrap replicates.

The script reads `DF/MIMIC/mimic31_aki_admission_base_v3.rds`, the landmark-aligned `LANDMARK_pred_dt.csv` predictions, and longitudinal SCr measurements under the project data directory. Analyses and tables are written to `Enrichment_analysis/` beneath that directory.

```bash
Rscript 06_Enrichment_analysis.R --base_dir=/secure/TrajAKI_Ventil
```
