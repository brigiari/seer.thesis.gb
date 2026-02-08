# ─────────────────────────────────────────────────────────────────────────────
# 04-analytic-dataset.R
# Select final model variables into `surv` dataset
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "03-variable-recoding.R"))

cat("\n── Building analytic dataset ─────────────────────────────────────\n")

surv <- nsclc_incl |>
  select(
    # treatment
    radiation_recode,
    chemotherapy_recode_yes_no_unk,
    surgery,
    reg_lymph_surg,

    # histology & tumor characteristics
    aya_site_recode_2020_revision,
    grade_recode_thru_2017,
    histology_recode_broad_groupings,
    visceral_and_parietal_pleural_invasion_recode_2010,
    separate_tumor_nodules_ipsilateral_lung_recode_2010,

    # demographics
    year_of_diagnosis,
    diagnostic_confirmation,
    age,
    sex,
    marital_status_at_diagnosis,
    household_income,
    rural_urban,
    race_recode_white_black_other,
    laterality,
    primary_site_labeled,

    # nodes
    regional_nodes_examined,
    regional_nodes_positive,

    # staging
    combined_summary_stage_2004,
    cs_site_specific_factor_1_2004_2017_varying_by_schema,
    cs_site_specific_factor_2_2004_2017_varying_by_schema,

    # distant metastases
    seer_combined_mets_at_dx_bone_2010,
    seer_combined_mets_at_dx_brain_2010,
    seer_combined_mets_at_dx_liver_2010,

    # AJCC 7th edition TNM
    derived_ajcc_t_7th_ed_2010_2015,
    derived_ajcc_n_7th_ed_2010_2015,
    derived_ajcc_m_7th_ed_2010_2015,

    # outcome
    survival_months,
    cs_status
  )

cat("  Analytic dataset `surv`:", format(nrow(surv), big.mark = ","),
    "records x", ncol(surv), "variables\n")

cat("\n── CONSORT flow (final) ──────────────────────────────────────────\n")
print(flow, n = Inf)
