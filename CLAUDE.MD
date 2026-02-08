# SEER Lung Cancer — Random Survival Forest Analysis

## Project Overview
Scientific manuscript analyzing NSCLC (non-small cell lung cancer) survival using SEER registry data and Random Survival Forests (RSF). The goal is to measure RSF performance and utility for cancer-specific survival prediction, with SHAP-based interpretability.

## Data
- **Source**: SEER Research Plus, 17 registries, Nov 2023 submission (2000–2021 diagnoses, follow-up through 2022)
- **Raw data path**: stored externally via `Sys.getenv("RAW_PATH")` (set in `.Renviron`); session export folder is `/data-raw/sess-10`
- **Size**: ~100k cases before filtering; final analytic cohort is substantially smaller after inclusion/exclusion criteria
- **Outcome**: Cancer-specific death (`cs_status`, binary) with `survival_months` as time variable

## Cohort Definition
1. **Site**: Lung and Bronchus (ICD-O-3 site recode)
2. **Primary site**: C34.0–C34.3, C34.9
3. **Behavior**: Malignant only (behavior code 3)
4. **Histology**: NSCLC codes (adenocarcinoma, squamous, large-cell, other specified); excludes small-cell and NOS
5. **Diagnostic confirmation**: Microscopically confirmed (values 1–4)
6. **Staging**: AJCC 7th edition; excludes T0, TX, T1NOS, T2NOS, NX, M1/M1NOS, OCCULT, UNK Stage
7. **Other**: Primary tumors only, age ≥ 18, excludes death-certificate-only and autopsy-only

## Key Variables
| Category | Variables |
|---|---|
| **Outcome** | `survival_months`, `cs_status` |
| **Treatment** | `surgery`, `radiation_recode`, `chemotherapy_recode_yes_no_unk`, `reg_lymph_surg` |
| **Staging** | `derived_ajcc_t/n/m_7th_ed_2010_2015`, `combined_summary_stage_2004` |
| **Tumor** | `histology_recode_broad_groupings`, `grade_recode_thru_2017`, `visceral_and_parietal_pleural_invasion_recode_2010`, `separate_tumor_nodules_ipsilateral_lung_recode_2010`, `primary_site_labeled` |
| **Demographics** | `age`, `sex`, `race_recode_white_black_other`, `marital_status_at_diagnosis` |
| **Socioeconomic** | `household_income`, `rural_urban` |
| **Nodes** | `regional_nodes_examined`, `regional_nodes_positive` |
| **Metastases** | `seer_combined_mets_at_dx_bone/brain/liver_2010` |
| **CS factors** | `cs_site_specific_factor_1/2_2004_2017_varying_by_schema` |

## Modeling Approach
- **Algorithm**: Random Survival Forest via `randomForestSRC`
- **Tuning**: `tune()` with 200 trees → optimal `mtry`, `nodesize`
- **Final model**: 1000 trees, OOB-based evaluation
- **Variable selection**: Permutation VIMP, threshold > 0.001
- **Missing data**: Internal forest-based imputation (`na.action = "na.impute"`)
- **Metrics**: C-index (Harrell), Brier score (KM and RSF censoring), CRPS
- **Interpretability**: SHAP values via `kernelshap` + `shapviz`

## Clinician Review Notes (MVV)
Key limitations to address in manuscript:
1. **Pleura invasion bias**: ~64% unknown — likely non-operable patients, unfavorable prognosis bias
2. **Missing anatomopathological variables**: LVI and perineural invasion unavailable in SEER
3. **Surgery-mortality bias**: ~63% mortality likely concentrated in non-operated patients
4. **Treatment sequencing**: Binary yes/no treatment variables oversimplify; no sequence variables available in this data extract

## Coding Conventions
- **Language**: R (≥ 4.4)
- **Style**: tidyverse preferred (pipes `|>` or `%>%`, dplyr verbs)
- **Reproducibility**: `set.seed(666)` throughout
- **Project management**: RStudio project (`seer.thesis.gb.Rproj`), `here::here()` for paths
- **Data privacy**: Raw data path in `.Renviron` (gitignored), no data committed to repo

## Directory Structure
```
seer.thesis.gb/
├── background/           # Clinician comments PDF, SEER resources, staging manuals
│   ├── MVV comments eng.pdf
│   └── Resources/        # SEER dictionaries, recoding guides, staging PDFs
├── script/
│   ├── final/            # Current working scripts (canonical)
│   └── ...               # Legacy/exploratory scripts
├── output/               # Reports, figures
├── .Renviron             # RAW_PATH env var (gitignored)
└── seer.thesis.gb.Rproj
```

## Important Notes
- Data is very large (~100k); use subsampling (~10k) for development, full dataset for final runs
- Age filter: 18–80 (resolved)
- Include 0-month survivors (resolved)
- Normalize/dummy-code only if needed — RSF handles factors natively
- `randomForestSRC` supports parallelization: `options(rf.cores = detectCores())`
- SHAP computation is expensive; sample ≤1000 rows for `kernelshap`
