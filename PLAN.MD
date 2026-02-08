# PLAN: SEER Lung Cancer RSF — Clean Rebuild

Rebuild the analysis from scratch with clean, modular scripts incorporating clinician (MVV) feedback, starting from the existing messy codebase.

## Phase 1: Project Setup

- [x] **1.1** Create clean directory structure: `script/new/`, `output/figures/`, `output/tables/`
- [x] **1.2** Write `script/new/00-setup.R` — load libraries, set seed, set parallelization, define paths
- [x] **1.3** Verify `.Renviron` has `RAW_PATH` set correctly

## Phase 2: Data Import & Cleaning

- [x] **2.1** Write `script/new/01-data-import.R` — load raw SEER data from `.RData`
- [x] **2.2** Write `script/new/02-cohort-selection.R` — apply all inclusion/exclusion filters:
  - Site recode: Lung and Bronchus
  - Primary site: C34.0–C34.3, C34.9
  - Behavior code: Malignant (3)
  - Histology: NSCLC codes (updated per MVV review)
  - Diagnostic confirmation: Microscopically confirmed (1–4)
  - Exclude NOS histologies
  - Primary tumors only, age ≥ 18, exclude death-cert/autopsy-only
  - AJCC 7th ed filters (exclude T0, TX, NX, M1, OCCULT, UNK)
  - Drop NAs in staging, survival, event status
- [x] **2.3** Write `script/new/03-variable-recoding.R` — recode and engineer features:
  - `cs_status` (binary event), `age` (numeric), `surgery`, `systemic_therapy`
  - `marital_status_at_diagnosis` (binary), `household_income` (grouped), `rural_urban`
  - `reg_lymph_surg` (grouped), `regional_nodes_examined/positive` (cleaned)
  - Drop zero-variance and all-NA columns
- [x] **2.4** Write `script/new/04-analytic-dataset.R` — select final model variables into `surv` dataset
- [x] **2.5** Add flowchart counts: log N at each filtering step for CONSORT-style diagram

## Phase 3: Descriptive Analysis

- [x] **3.1** Write `script/new/05-descriptive.R` — `gtsummary::tbl_summary()` for full cohort
- [x] **3.2** Create surgery-mortality cross-table (MVV comment: 63% mortality bias)
- [x] **3.3** Tabulate missingness in pleura invasion variable (MVV comment: 64% unknown)
- [ ] **3.4** Explore treatment sequence variables available in SEER (document for limitations)

## Phase 4: Modeling — Development Sample (~10k)

- [x] **4.1** Write `script/new/06-sampling.R` — stratified subsample ~10k for development
- [x] **4.2** Write `script/new/07-rsf-tune.R` — hyperparameter tuning (`tune()`, 200 trees)
- [x] **4.3** Write `script/new/08-rsf-train.R` — train RSF (1000 trees, tuned params, OOB)
- [x] **4.4** Permutation VIMP + variable selection (threshold > 0.001) — in `script/new/08-rsf-train.R`
- [x] **4.5** Refit reduced model with selected variables — in `script/new/08-rsf-train.R`
- [x] **4.6** Write `script/new/09-rsf-metrics.R` — C-index, Brier score, CRPS, error rate plots

## Phase 5: Modeling — Full Dataset

- [ ] **5.1** Re-run tuning + training on full dataset with same pipeline
- [ ] **5.2** Complete-case sensitivity analysis (train/test 70/30 split)
- [ ] **5.3** Compare splitting rules (logrank, bs.gradient, logrankscore)
- [ ] **5.4** Nested variable importance analysis (incremental error reduction)

## Phase 6: Interpretability

- [ ] **6.1** Write `script/new/10-shap.R` — SHAP values via `kernelshap` (sample ≤1000 rows)
- [ ] **6.2** Generate SHAP plots: beeswarm importance, bar importance, waterfall, force
- [ ] **6.3** Marginal and partial dependence plots for top predictors
- [ ] **6.4** Minimal depth variable selection as supplementary analysis

## Phase 7: Figures & Tables for Manuscript

- [ ] **7.1** Publication-quality survival curves (OOB)
- [ ] **7.2** VIMP plot with confidence intervals (subsample/jackknife)
- [ ] **7.3** Brier score and CRPS plots
- [ ] **7.4** SHAP summary figures
- [ ] **7.5** Descriptive Table 1
- [ ] **7.6** Surgery-mortality bias table (MVV request)
- [ ] **7.7** Save all figures to `output/figures/` at 300 dpi

## Phase 8: Manuscript Support

- [ ] **8.1** Write limitations paragraph addressing:
  - Pleura invasion bias (64% unknown)
  - Missing LVI and perineural invasion
  - Surgery-mortality confounding (63% mortality)
  - Treatment sequencing oversimplification (binary yes/no)
- [ ] **8.2** Update `output/report-rsf.Rmd` to source new clean scripts
- [ ] **8.3** Final reproducibility check: clean R session, run all scripts end-to-end

## Decisions (Resolved)

1. ~~Include treatment sequence variables?~~ → **No**, not available in data. Ignore.
2. ~~Include pleura invasion variable despite 64% missingness?~~ → **Yes**, keep for now.
3. ~~Age filtering?~~ → **18–80**.
4. ~~Survival time floor?~~ → **Include 0-month** survivors.
5. ~~Normalize / dummy-code?~~ → **Only if needed** (RSF handles factors natively).
