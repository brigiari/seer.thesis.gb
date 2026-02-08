library(tidyverse)
library(stringr)
library(gtsummary)
folder = "/data-raw/sess-10"  

load(file = paste0(Sys.getenv("RAW_PATH"), folder, "/data.RData"))

# define the NSCLC codes

nsclc_codes <- sort(c(
  # Large-cell carcinomas
  8012:8014, 8021, 8034, 8082,              
  # Squamous cell carcinomas
  8051:8052, 8070:8076, 8078, 8083:8084, 8090, 8094, 8123 ,              
  # Adenocarcinomas
  8015, 8050, 8140:8141, 8143:8145, 8147, 8190, 8201, 8211, 
  8250:8255, 8260, 8290, 8310, 8320, 8323, 8333, 8401, 8440, 
  8470:8471, 
  8480:8481, 8490, 8503, 8507, 8550, 8570:8572, 8574, 8576            
  # NOS
  # 8046, 8003:8004, 8022, 8030, 8031:8033, 8035, 8120, 8200, 
  # 8240:8241, 8243:8246, 8249, 8430, 8525, 8560, 8562, 8575                   
))

small_cell_codes <- c(8002, 8041:8045)

# 1. Initial filtering & basic recoding -----------------------------------------

db <- data |>
  # filter only primary tumor
  filter(primary_by_international_rules == "Yes") |>
  
  # Keep only lung & bronchus cases (double checked also with primary_site_labeled)
  filter(site_recode_icd_o_3_2023_revision == "Lung And Bronchus") |>  
  
  filter(
    histologic_type_icd_o_3 %in% nsclc_codes # nsclc_codes, nsclc_codes
               ) |>  
  
  
  # Remove any records missing T, N, M stage
  filter(!is.na(derived_ajcc_t_7th_ed_2010_2015)) |>  
  filter(!is.na(derived_ajcc_n_7th_ed_2010_2015)) |>   
  filter(!is.na(derived_ajcc_m_7th_ed_2010_2015)) |> 
  filter(!is.na(derived_ajcc_stage_group_7th_ed_2010_2015)) |> 
  
  
  
  # Convert many character columns to factors, recode numeric date- and size-fields
  mutate(
    across(where(is.character), factor),                                    # all character → factor
    time_from_diagnosis_to_treatment_in_days_recode = 
      as.numeric(time_from_diagnosis_to_treatment_in_days_recode),         # days → numeric
    survival_months = as.numeric(survival_months),                         # survival time → numeric
    total_number_of_in_situ_malignant_tumors_for_patient = 
      as.numeric(total_number_of_in_situ_malignant_tumors_for_patient),    # count → numeric
    tumor_size_over_time_recode_1988 = 
      as.numeric(tumor_size_over_time_recode_1988),                       # size → numeric
    # create binary event indicator: 1 = cancer-specific death, 0 = other/ alive
    cs_status = if_else(
      seer_cause_specific_death_classification == 
        "Dead (attributable to this cancer dx)", 1, 0
    ),
    # strip “ years” and convert age to numeric
    age = as.numeric(str_remove_all(age_recode_with_single_ages_and_90, " years"))
  ) |> 
  
  filter(age %in% c(18:80)) %>% 
  
  filter(!type_of_reporting_source %in% c(
    "Death certificate only",
    "Autopsy only"
  )) |> 

  # ───────────────────────────────────────────────────────────────────────────
  # 2. Drop unwanted columns & rows with all NAs
  # ───────────────────────────────────────────────────────────────────────────

  select(
    -where(~ all(is.na(.))),               # drop columns entirely NA
    -primary_site,                         # drop redundant site codes
    # -histologic_type_icd_o_3, 
    # -icd_o_3_hist_behav,
    -survival_months_flag, 
    -patient_id,                           # drop identifier
    # -age_recode_with_1_year_olds,          # drop alternate age fields
    # -age_recode_with_1_year_olds_and_90,
    # -age_recode_with_single_ages_and_85,
    # -tumor_size_over_time_recode_1988,     # original size variable now numeric
    # -contains("cod"),                      # any column with “cod” in its name
    -seer_registry_with_ca_and_ga_as_whole_states
  ) |> 
  drop_na(survival_months, cs_status) |>     # remove any remaining NAs in time/event
  filter(survival_months != 0) 



# ─────────────────────────────────────────────────────────────────────────────
# 3. Quick table of all variables (counts & missingness)
# ─────────────────────────────────────────────────────────────────────────────

# db |>
#   gtsummary::tbl_summary(
#     missing = "ifany"   # show n (%) missing if any
#   ) |>
#   gtsummary::add_n()    # add total N column


# ─────────────────────────────────────────────────────────────────────────────
# 4. Further recoding for survival modeling
# ─────────────────────────────────────────────────────────────────────────────

surv <- db |> 
  filter(diagnostic_confirmation%in%c(
    "Positive histology",
    "Pos hist AND immunophenotyping AND/OR pos genetic studies")
    ) |> 
  mutate(
    # simplify marital status: married → “1”, all else → “0”
    marital_status_at_diagnosis = if_else(
      marital_status_at_diagnosis == "Married (including common law)", "1", "0"
    ),
    
    # collapse diagnostic confirmation into Clinical vs Pathological vs Unknown
    # diagnostic_confirmation = case_when(
    #   diagnostic_confirmation %in% c(
    #     "Clinical diagnosis only",
    #     "Direct visualization without microscopic confirmation",
    #     "Clinical diagnosis only, direct visualization without microscopic confirmation",
    #     "Radiography without microscopic confirm"
    #   ) ~ "Clinical",
    #   diagnostic_confirmation == "Unknown" ~ "Unknown",
    #   TRUE ~ "Pathological"
    # ),
    
    # unify laterality: any indication of side → “1”, else → “2”
    # laterality = case_when(
    #   laterality %in% c(
    #     "Left - origin of primary",
    #     "Right - origin of primary",
    #     "Only one side - side unspecified",
    #     "Not a paired site"
    #   ) ~ "1",
    #   TRUE ~ "2"
    # ),
    
    # extract numeric T-category (T1–T4) or keep as is
    cT = case_when(
      str_detect(derived_ajcc_t_7th_ed_2010_2015, "T1") ~ "T1",
      str_detect(derived_ajcc_t_7th_ed_2010_2015, "T2") ~ "T2",
      str_detect(derived_ajcc_t_7th_ed_2010_2015, "T3") ~ "T3",
      str_detect(derived_ajcc_t_7th_ed_2010_2015, "T4") ~ "T4",
      TRUE ~ derived_ajcc_t_7th_ed_2010_2015
    ),
    
    # N and M are taken directly from AJCC recode
    cN = derived_ajcc_n_7th_ed_2010_2015,
    M  = derived_ajcc_m_7th_ed_2010_2015,
    
    # collapse detailed stage into I, II, III, IV using regex
    stage = case_when(
      # I or I + letter that is not I/V
      str_detect(
        derived_ajcc_stage_group_7th_ed_2010_2015,
        regex("^I(?:$|[^IV][A-Za-z]*)$")
      ) ~ "I",
      # II but not III, IIV, etc.
      str_detect(
        derived_ajcc_stage_group_7th_ed_2010_2015,
        regex("^II(?!I)[A-Za-z]*$")
      ) ~ "II",
      # III but not IIII
      str_detect(
        derived_ajcc_stage_group_7th_ed_2010_2015,
        regex("^III(?!I)[A-Za-z]*$")
      ) ~ "III",
      # IV but not IVI, IVV
      str_detect(
        derived_ajcc_stage_group_7th_ed_2010_2015,
        regex("^IV(?!I)[A-Za-z]*$")
      ) ~ "IV",
      # fallback: leave original
      TRUE ~ derived_ajcc_stage_group_7th_ed_2010_2015
    ),
    
    
    stage = if_else(
      stage %in% c("OCCULT", "UNK Stage"), NA_character_, stage
    ),  # convert "Unknown" to NA
    
    cT = if_else(cT%in%c("T0", "TX"), NA_character_, cT),
    cN = if_else(cN=="NX", NA_character_, cN),
    
    
    
    regional_nodes_positive_1988 = case_when(
      regional_nodes_positive_1988 %in% c(1:90) ~ as.numeric(regional_nodes_positive_1988),
      regional_nodes_positive_1988 == 0 ~ 0,
      regional_nodes_positive_1988 %in% c(95, 97) ~ 1,
      regional_nodes_positive_1988 %in% c(98, 99) ~ NA_real_,
      TRUE ~ regional_nodes_positive_1988
    ),
    
    regional_nodes_positive = case_when(
      regional_nodes_positive_1988 <= 20 ~ "<20",
      regional_nodes_positive_1988 <= 40 ~ "21-40",
      regional_nodes_positive_1988 <= 60 ~ "41-60",
      TRUE ~ "60+"
    ),
    
    regional_nodes_examined_1988 = case_when(
      regional_nodes_examined_1988 %in% c(1:90) ~ as.numeric(regional_nodes_examined_1988),
      regional_nodes_examined_1988 == 0 ~ 0,
      regional_nodes_examined_1988 %in% c(95) ~ 0,
      regional_nodes_examined_1988 %in% c(96:99) ~ NA_real_,
      TRUE ~ regional_nodes_examined_1988
    ),
    
    regional_nodes_examined = case_when(
      regional_nodes_examined_1988 <= 20 ~ "<20",
      regional_nodes_examined_1988 <= 40 ~ "21-40",
      regional_nodes_examined_1988 <= 60 ~ "41-60",
      TRUE ~ "60+"
    ),
    
    surgery = case_when(
      rx_summ_surg_prim_site_1998 == "0" ~ "None",
      rx_summ_surg_prim_site_1998 %in% c(10:19) ~ "Tumor destruction",
      rx_summ_surg_prim_site_1998 %in% c(20:80) ~ "Resection",
      rx_summ_surg_prim_site_1998 == 90 ~ "Surgery performed, type unknown",
      rx_summ_surg_prim_site_1998 == "99" ~ "Unknown",
      
      TRUE ~ "Unknown"
    ),
    
    reg_lymph_surg = case_when(
      rx_summ_scope_reg_ln_sur_2003 == "None" ~ "None removed",
      rx_summ_scope_reg_ln_sur_2003 == "Sentinel lymph node biopsy" ~ "Only sentinel removed",
      rx_summ_scope_reg_ln_sur_2003 %in% c("1 to 3 regional lymph nodes removed",
                                           "Biopsy or aspiration of regional lymph node, NOS",
                                           "4 or more regional lymph nodes removed",
                                           "Number of regional lymph nodes removed unknown") ~ "Regional lymph nodes removed",
      rx_summ_scope_reg_ln_sur_2003 %in% c("Sentinel node biopsy and lym nd removed different times",
                                           "Sentinel node biopsy and lym nd removed same/unstated time") ~ "Regional + Sentinel",
      rx_summ_scope_reg_ln_sur_2003 == "Unknown or not applicable" ~ "Unknown",
      TRUE ~ "Unknown"
    ),
    
    radiation = if_else(radiation_recode %in% 
                          c("None/Unknown",
                            "Recommended, unknown if administered",
                            "Refused (1988+)"), "No/Unknown", "Yes"),
    
    age = case_when(
      age <= 40 ~ "18-40",
      age <= 60 ~ "41-60",
      TRUE ~ "61+"
    ),
    
    household_income = case_when(
      median_household_income_inflation_adj_to_2022 %in% c(
        "< $40,000", "$40,000 - $44,999", "$45,000 - $49,999", 
        "$50,000 - $54,999",  "$55,000 - $59,999"
      ) ~ "< $60,000",
      
      median_household_income_inflation_adj_to_2022 %in% c(
        "$60,000 - $64,999", "$65,000 - $69,999", "$70,000 - $74,999",
        "$75,000 - $79,999", "$80,000 - $84,999") ~ "$60,000 - $84,999",
      
      median_household_income_inflation_adj_to_2022 == "$120,000+" ~ "$120,000+",
      TRUE ~ "$85,000 - $119,099"),
    
    rural_urban = case_when(
      str_detect(rural_urban_continuum_code, "Nonmetropolitan") ~ "Non-metropolitan",
      str_detect(rural_urban_continuum_code, " metropolitan") ~ "Metropolitan",
      TRUE ~ NA_character_
    ),
        
    
    
    
    # radiation = if_else(rx_summ_surg_rad_seq == "No radiation and/or no surgery; unknown if surgery and/or radiation given",
    #                     "None", "Yes"),
    
    # systemic_th = if_else(rx_summ_systemic_sur_seq_2007 == "No systemic therapy and/or surgical procedures",
    #                       "None", "Yes"),
    
    # re-factor any remaining character cols
    across(where(is.character), factor)
    
    
  ) |> 
  
  # ───────────────────────────────────────────────────────────────────────────
  # 5. Final cleanup & select modeling variables
  # ───────────────────────────────────────────────────────────────────────────
  
  drop_na(cT, cN, stage) |>       # drop T0 if present
  
  select(
    radiation,
    chemotherapy_recode_yes_no_unk,
    # histology_recode_broad_groupings, # add for NSCLC
    aya_site_recode_2020_revision,
    grade_recode_thru_2017,
    # separate_tumor_nodules_ipsilateral_lung_recode_2010,
    visceral_and_parietal_pleural_invasion_recode_2010,
    # year_of_diagnosis,


    # core demographics & diagnosis
    age, 
    sex,
    primary_site_labeled,
    marital_status_at_diagnosis,
    household_income,
    rural_urban,
    race_recode_white_black_other,


    # icd_o_3_hist_behav_malignant,
    primary_site_labeled,
    # treatment variables
    surgery, reg_lymph_surg, 
    # radiation, systemic_th,
    # cs variables
    # cs_lymph_nodes_2004_2015:cs_lymph_nodes_2004_2015,
    # cs_tumor_size_ext_eval_2004_2015:cs_mets_eval_2004_2015,
    
    # regional_nodes_positive,
    regional_nodes_examined,
    regional_nodes_positive,
    # cs_site_specific_factor_1_2004_2017_varying_by_schema,
    # cs_site_specific_factor_2_2004_2017_varying_by_schema,


    # outcome
    # seer_cause_specific_death_classification,
    survival_months, cs_status,

    # AJCC-derived TNM & stage
    cT, cN, M, stage
  )  

