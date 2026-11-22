# Predicting Heart Failure Mortality Using Clinical Biomarkers

**Author:** John Ryan Whetzel  
**Course:** DATA 230 — Phase 4 Reproducible Project Package

## Overview

This project builds and evaluates a binomial logistic regression model to predict mortality risk (`DEATH_EVENT`) among heart failure patients using the **Heart Failure Clinical Records Dataset** ($N = 299$, UCI Machine Learning Repository, ID: 519). The workflow uses a stratified 70/30 train/test split (`seed = 230`), a `0.50` classification probability threshold, five clinical predictors (`age`, `ejection_fraction`, `serum_creatinine`, `serum_sodium`, and `time`), a training-prevalence baseline comparison, and Cook's distance influence diagnostics.

## How to Run

1. Open a fresh R session and set the working directory to the project root folder.
2. Ensure internet access to download the dataset from UCI, or place `heart_failure_clinical_records_dataset.csv` inside the project root directory.
3. Run the master script from the R console:
   ```r
   source("final_analysis.R")
   ```
   (You can also run `Rscript final_analysis.R` from a terminal in the project root.)
4. Inspect the generated metrics, confusion matrix, figures, and software version log in their respective output folders.

## Execution Order

This workflow runs through a single self-contained script:
1. `final_analysis.R` — Checks for and installs required packages if missing, loads or downloads the dataset, formats `DEATH_EVENT` (`survived` vs. `deceased`), creates the stratified 70/30 split (`seed = 230`), fits the logistic regression model, computes held-out test metrics and baseline comparisons, generates diagnostic and performance plots, and writes R and package versions to disk.

## Required Packages

- `readr`
- `dplyr`
- `tibble`
- `ggplot2`
- `pROC`

## Software Versions

Exact versions of R and all loaded namespaces from your run are written automatically to:
- `phase4_work/analysis_software_versions.txt`

## Data Access

- **Dataset:** Heart Failure Clinical Records Dataset (UCI Machine Learning Repository, ID: 519)
- **Source URL:** https://archive.ics.uci.edu/dataset/519/heart+failure+clinical+records
- **Local Fallback:** If automatic download is unavailable, place `heart_failure_clinical_records_dataset.csv` anywhere inside the project directory and run `final_analysis.R` again.

## Recreated Outputs

Running `final_analysis.R` in a clean session recreates:
- `phase2_work/cleaned_data.csv` (formatted dataset cache)
- `phase3_work/metrics_table.csv` (held-out accuracy, precision, recall, and ROC-AUC for the logistic model vs. training-prevalence baseline)
- `phase3_work/confusion_matrix.csv` (held-out confusion matrix at threshold `0.50`)
- `phase3_work/figs/fig_01_roc.png` (held-out ROC curve)
- `phase3_work/figs/fig_02_confusion_matrix.png` (held-out confusion matrix plot)
- `phase3_work/figs/fig_03_training_influence.png` (Cook's distance training influence plot with `4/n` reference threshold)
- `phase4_work/analysis_software_versions.txt` (R and package version log)

## Special Notes & Limitations

- Always run from the project root directory in a fresh session so relative paths (`phase2_work/`, `phase3_work/`, `phase4_work/`) resolve cleanly.
- Training records above the `4/n` Cook's distance threshold reflect clinically severe profiles and are kept in the model rather than dropped automatically.
- With a moderate sample size ($N = 299$) and baseline clinical variables, a `0.50` cutoff still produces false negatives, so predictions should support—not replace—clinical judgment.

## References and Tools

- Chicco, D., & Jurman, G. (2020). Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. *BMC Medical Informatics and Decision Making*, 20(1), 16. https://doi.org/10.1186/s12911-020-1023-5
- UCI Machine Learning Repository. (2020). *Heart failure clinical records* [Dataset]. https://archive.ics.uci.edu/dataset/519/heart+failure+clinical+records
- R Core Team. (2026). *R: A language and environment for statistical computing* [Computer software]. R Foundation for Statistical Computing. https://www.R-project.org/
- Wickham, H., et al. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686. https://doi.org/10.21105/joss.01686
- Robin, X., et al. (2011). pROC: An open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*, 12, 77. https://doi.org/10.1186/1471-2105-12-77
- Google. (2026). *Gemini* [Large language model]. Used for workflow structuring, code debugging, and drafting support.