# ─────────────────────────────────────────────────────────────────────────────
# 03-variable-recoding.R
# Recode and engineer features from the filtered cohort
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "02-cohort-selection.R"))

cat("\n── Variable recoding ─────────────────────────────────────────────\n")

nsclc_incl <- cohort |>
  mutate(
    # marital status → binary
    marital_status_at_diagnosis = if_else(
      marital_status_at_diagnosis == "Married (including common law)", "1", "0"
    ),

    # regional nodes positive — clean sentinel/unknown codes
    regional_nodes_positive = case_when(
      regional_nodes_positive_1988 %in% c(95, 97, 98, 99, 126) ~ NA_real_,
      TRUE ~ regional_nodes_positive_1988
    ),

    # regional nodes examined — clean codes ≥ 95
    regional_nodes_examined = case_when(
      regional_nodes_examined_1988 %in% c(0:90) ~ as.numeric(regional_nodes_examined_1988),
      regional_nodes_examined_1988 >= 95 ~ NA_real_
    ),

    # surgery — from reason_no_cancer_directed_surgery
    surgery = case_when(
      reason_no_cancer_directed_surgery == "Surgery performed" ~ "Yes",
      reason_no_cancer_directed_surgery %in% c(
        "Unknown; death certificate; or autopsy only (2003+)",
        "Recommended, unknown if performed"
      ) ~ NA_character_,
      TRUE ~ "No"
    ),

    # systemic therapy — from rx_summ_systemic_sur_seq_2007
    systemic_therapy = case_when(
      rx_summ_systemic_sur_seq_2007 %in%
        c("No systemic therapy and/or surgical procedures") ~ "No",
      rx_summ_systemic_sur_seq_2007 %in% c("Sequence unknown") ~ NA_character_,
      TRUE ~ "Yes"
    )
  ) |>

  drop_na(regional_nodes_examined_1988, surgery, systemic_therapy) |>
  log_step("Drop NA in nodes_examined / surgery / systemic_therapy")

nsclc_incl <- nsclc_incl |>
  mutate(
    # regional lymph node surgery — grouped
    reg_lymph_surg = case_when(
      rx_summ_scope_reg_ln_sur_2003 == "None" ~ "None removed",
      rx_summ_scope_reg_ln_sur_2003 == "Sentinel lymph node biopsy" ~ "Only sentinel removed",
      rx_summ_scope_reg_ln_sur_2003 %in% c(
        "1 to 3 regional lymph nodes removed",
        "Biopsy or aspiration of regional lymph node, NOS",
        "4 or more regional lymph nodes removed",
        "Number of regional lymph nodes removed unknown"
      ) ~ "Regional lymph nodes removed",
      rx_summ_scope_reg_ln_sur_2003 %in% c(
        "Sentinel node biopsy and lym nd removed different times",
        "Sentinel node biopsy and lym nd removed same/unstated time"
      ) ~ "Regional + Sentinel",
      rx_summ_scope_reg_ln_sur_2003 == "Unknown or not applicable" ~ "Unknown",
      TRUE ~ "Unknown"
    ),

    # household income — grouped into 3 bands
    household_income = case_when(
      median_household_income_inflation_adj_to_2022 %in% c(
        "< $40,000", "$40,000 - $44,999", "$45,000 - $49,999",
        "$50,000 - $54,999", "$55,000 - $59,999"
      ) ~ "< $60,000",
      median_household_income_inflation_adj_to_2022 %in% c(
        "$60,000 - $64,999", "$65,000 - $69,999", "$70,000 - $74,999",
        "$75,000 - $79,999", "$80,000 - $84,999"
      ) ~ "$60,000 - $84,999",
      median_household_income_inflation_adj_to_2022 == "$120,000+" ~ "$120,000+",
      TRUE ~ "$85,000 - $119,099"
    ),

    # rural/urban — simplified
    rural_urban = case_when(
      str_detect(rural_urban_continuum_code, "Nonmetropolitan") ~ "Non-metropolitan",
      str_detect(rural_urban_continuum_code, " metropolitan") ~ "Metropolitan",
      TRUE ~ NA_character_
    ),

    # re-factor any remaining character cols and drop unused levels
    across(where(is.character), factor),
    across(where(is.factor), ~ fct_drop(.))
  ) |>

  # remove all-NA columns
  select(-where(~ all(is.na(.)))) |>
  log_step("After variable recoding + drop all-NA columns")

cat("  Final recoded cohort:", format(nrow(nsclc_incl), big.mark = ","), "records\n")
