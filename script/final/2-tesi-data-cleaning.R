library(tidyverse)
library(stringr)
library(gtsummary)
folder = "/data-raw/sess-10"  
# usethis::edit_r_environ("project")
load(file = paste0(Sys.getenv("RAW_PATH"), folder, "/datav2.RData"))

# define the NSCLC codes

nsclc_codes <- sort(c(
  # Large-cell carcinomas
  8012:8014, 8021, 8034, 8082,              
  # Squamous cell carcinomas
  8051:8052, 8070:8076, 8078, 8083:8084, 8090, 8094, 8123 ,              
  # Adenocarcinomas
  8015, 8050, 8140:8141, 8143:8145, 8147, 8190, 8201, 8211, 
  8250:8255, 8260, 8290, 8310, 8323, 8333, 8401, 8440, 
  8470:8471, 
  8480:8481, 8490, 8503, 8507, 8550, 8570:8572, 8574, 8576,            
  # NOS
  8046, 8003:8004, 8022, 8030, 8031:8033, 8035, 8120, 8200,
  8240:8241, 8243:8246, 8249, 8430, 8525, 8560, 8562, 8575
))

small_cell_codes <- c(8002, 8041:8045)
nos <- c(8046, 8003:8004, 8022, 8030, 8031:8033, 8035, 8120, 8200,
         8240:8241, 8243:8246, 8249, 8430, 8525, 8560, 8562, 8575)


nsclc <- data |> 
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
  
  filter(site_recode_icd_o_3_2023_revision == "Lung And Bronchus") |>
  filter(str_detect(primary_site_labeled, "34.0|34.1|34.2|34.3|34.8|34.9")) |> 
  
  filter(
    histologic_type_icd_o_3 %in% nsclc_codes # nsclc_codes, nsclc_codes
  )

#NSCLC patients in SEER database
nrow(nsclc)
summary(nsclc$year_of_diagnosis)

table(nsclc$primary_site_labeled, nsclc$laterality)
table(nsclc$rx_summ_surg_prim_site_1998)

nsclc |> 
  filter(
    histologic_type_icd_o_3 %in% nos # nsclc_codes, nsclc_codes
  ) |> nrow()

nsclc_incl <- nsclc |> 
  filter(
    !histologic_type_icd_o_3 %in% nos # nsclc_codes, nsclc_codes
  ) |> 
  
  # ───────────────────────────────────────────────────────────────────────────
  # Some filtering
  # ───────────────────────────────────────────────────────────────────────────
  
  # only ones with + histology (as per IASLC suggestion)
  filter(diagnostic_confirmation %in% "Positive histology") |>
  
  # filter only primary tumor
  filter(primary_by_international_rules == "Yes") |>
  
  
  # filter by reporting source
  filter(!type_of_reporting_source %in% c(
    "Death certificate only",
    "Autopsy only"
  )) |>
  
  # filter by age
  filter(age >= 18) |>
  
  # # escludo i NOS, T0 e TX
  filter(!derived_ajcc_t_7th_ed_2010_2015 %in% c(
    "T0", "T1NOS", "T2", "T2NOS", "TX"
  )) |> 
  
  
  filter(!derived_ajcc_n_7th_ed_2010_2015 %in% c(
    "NX"
  )) |>
  
  filter(!derived_ajcc_m_7th_ed_2010_2015 %in% c(
    "M1NOS", "M1"
  )) |> 
  
  filter(!derived_ajcc_stage_group_7th_ed_2010_2015 %in% c(
    "OCCULT", "UNK Stage"
  )) |> 
  
  drop_na(
    derived_ajcc_t_7th_ed_2010_2015,
    derived_ajcc_n_7th_ed_2010_2015,
    derived_ajcc_m_7th_ed_2010_2015,
    derived_ajcc_stage_group_7th_ed_2010_2015
  ) |>
  
  mutate(
    marital_status_at_diagnosis = if_else(
      marital_status_at_diagnosis == "Married (including common law)", "1", "0"
    ),
    
    regional_nodes_positive_198 = case_when(
      regional_nodes_positive_1988 %in% c(1:90) ~ as.numeric(regional_nodes_positive_1988),
      regional_nodes_positive_1988 == 0 ~ 0,
      regional_nodes_positive_1988 %in% c(95, 97) ~ 1,
      regional_nodes_positive_1988 %in% c(98, 99) ~ NA_real_,
      TRUE ~ regional_nodes_positive_1988
    ),
    
    surgery = case_when(
      reason_no_cancer_directed_surgery == "Surgery performed" ~ "Yes",
      reason_no_cancer_directed_surgery %in% c(
        "Unknown; death certificate; or autopsy only (2003+)",
        "Recommended, unknown if performed") ~ NA_character_,
      TRUE ~ "No"
    ),
    
    systemic_therapy = case_when(
      rx_summ_systemic_sur_seq_2007 %in% c("No systemic therapy and/or surgical procedures") ~ "No",
      rx_summ_systemic_sur_seq_2007 %in% c("Sequence unknown") ~ NA_character_,
      TRUE ~ "Yes"
    ),
    
    regional_nodes_positive = case_when(
      regional_nodes_positive_1988 %in% c(95, 97, 98, 99, 126) ~ NA_real_,
      TRUE ~ regional_nodes_positive_1988
    ),
    
    regional_nodes_examined = case_when(
      regional_nodes_examined_1988 %in% c(0:90) ~ as.numeric(regional_nodes_examined_1988),
      regional_nodes_examined_1988 >= 95 ~ NA_real_
      )) |> 
  
  drop_na(regional_nodes_examined_1988, surgery, systemic_therapy) |>
  
  mutate(
  #   regional_nodes_examined = case_when(
  #     regional_nodes_examined_1988 <= 20 ~ "<20",
  #     regional_nodes_examined_1988 <= 40 ~ "21-40",
  #     regional_nodes_examined_1988 <= 60 ~ "41-60",
  #     TRUE ~ "60+"
  #   ),
    
    
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
    
    # radiation = if_else(rx_summ_surg_rad_seq %in% 
    #                       c("No radiation and/or no surgery; unknown if surgery and/or radiation given"),
    #                     "No", "Yes"),
    
    
    # age = case_when(
    #   age <= 40 ~ "18-40",
    #   age <= 60 ~ "41-60",
    #   TRUE ~ "61+"
    # ),
    
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
    
    # re-factor any remaining character cols
    across(where(is.character), factor),
    across(where(is.factor), ~fct_drop(.))) |> 
  
  
  drop_na(survival_months, cs_status) |>     # remove any remaining NAs in time/event
  # filter(survival_months >= 3) |> 
  select(-where(~all(is.na(.))))




nsclc_incl |>
  gtsummary::tbl_summary()
skimr::skim(nsclc_incl)





# create db for survival analysis
surv <- nsclc_incl |> 
  select(
    # radiation,
    radiation_recode,
    chemotherapy_recode_yes_no_unk,
    surgery, 
    # systemic_therapy,
    reg_lymph_surg, 
    
    
    # histology_recode_broad_groupings, # add for NSCLC
    aya_site_recode_2020_revision,
    grade_recode_thru_2017,
    histology_recode_broad_groupings,
    # separate_tumor_nodules_ipsilateral_lung_recode_2010,
    visceral_and_parietal_pleural_invasion_recode_2010,
    
    
    
    # core demographics & diagnosis
    year_of_diagnosis,
    diagnostic_confirmation,
    age_recode_with_1_year_olds, 
    sex,
    marital_status_at_diagnosis,
    household_income,
    rural_urban,
    race_recode_white_black_other,
    laterality,
    # icd_o_3_hist_behav_malignant,
    primary_site_labeled,
    # treatment variables
    # radiation, systemic_th,
    # cs variables
    # cs_lymph_nodes_2004_2015:cs_lymph_nodes_2004_2015,
    # cs_tumor_size_ext_eval_2004_2015:cs_mets_eval_2004_2015,
    
    # regional_nodes_positive,
    regional_nodes_examined,
    regional_nodes_positive,
    # rx_summ_scope_reg_ln_sur_2003,
    # rx_summ_surg_oth_reg_dis_2003,
    # rx_summ_surg_rad_seq,
    separate_tumor_nodules_ipsilateral_lung_recode_2010,
    # cs_site_specific_factor_1_2004_2017_varying_by_schema,
    # cs_site_specific_factor_2_2004_2017_varying_by_schema,
    combined_summary_stage_2004,
    # cs_lymph_nodes_2004_2015,
    # cs_mets_at_dx_2004_2015,
    # cs_mets_eval_2004_2015,
    # cs_reg_node_eval_2004_2015,
    # cs_tumor_size_2004_2015,
    # cs_tumor_size_ext_eval_2004_2015,
    cs_site_specific_factor_1_2004_2017_varying_by_schema,
    cs_site_specific_factor_2_2004_2017_varying_by_schema,
    
    
    seer_combined_mets_at_dx_bone_2010,
    seer_combined_mets_at_dx_brain_2010,
    seer_combined_mets_at_dx_liver_2010, #att! N/A or Unknown
    
    # outcome
    # seer_cause_specific_death_classification,
    survival_months, cs_status,
    
    # AJCC-derived TNM & stage
    derived_ajcc_t_7th_ed_2010_2015,
    derived_ajcc_n_7th_ed_2010_2015,
    derived_ajcc_m_7th_ed_2010_2015
    # derived_ajcc_stage_group_7th_ed_2010_2015
  )  


surv |>
  gtsummary::tbl_summary()

