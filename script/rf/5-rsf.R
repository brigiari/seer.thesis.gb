library(tidyverse)
folder = "/data-raw/sess-10"  

load(file = paste0(Sys.getenv("RAW_PATH"), folder, "/data.RData"))


db <- data |> 
  # filtro solo Lung
  filter(site_recode_icd_o_3_2023_revision == "Lung And Bronchus") |>

  filter(!is.na(derived_eod_2018_t_2018)) |>    
  mutate(
    across(where(is.character), factor),
    time_from_diagnosis_to_treatment_in_days_recode = as.numeric(time_from_diagnosis_to_treatment_in_days_recode),
    survival_months = as.numeric(survival_months),
    total_number_of_in_situ_malignant_tumors_for_patient = as.numeric(total_number_of_in_situ_malignant_tumors_for_patient),
    tumor_size_over_time_recode_1988 = as.numeric(tumor_size_over_time_recode_1988),
    cs_status = if_else(seer_cause_specific_death_classification == "Dead (attributable to this cancer dx)", 1, 0),
    age_recode_with_single_ages_and_90 = as.numeric(str_remove_all(age_recode_with_single_ages_and_90, " years")),
  ) |> 
  select(-where(~all(is.na(.))),
         -survival_months_flag, -patient_id,
         -age_recode_with_1_year_olds,
         -age_recode_with_1_year_olds_and_90,
         -age_recode_with_single_ages_and_85,
         -tumor_size_over_time_recode_1988,
         -contains("cod"),
         -seer_registry_with_ca_and_ga_as_whole_states) |> 
  drop_na(survival_months, cs_status) |> 
  filter(!survival_months == 0) 
  # select(derived_eod_2018_t_2018, 
  #        where(~ !(is.factor(.) && n_distinct(.) >= 10))) |> 
  # select(where(~ !(is.factor(.) && n_distinct(.) <= 1)))


db |> names()
surv <- db |> 
  
  mutate(
    marital_status_at_diagnosis = if_else(
      marital_status_at_diagnosis == "Married", 1, 0),
    
    diagnostic_confirmation = case_when(
      diagnostic_confirmation == "Clinical diagnosis only"~ "Clinical",
      diagnostic_confirmation == "Direct visualization without microscopic confirmation"~ "Clinical",
      diagnostic_confirmation == "Clinical diagnosis only, direct visualization without microscopic confirmation"~ "Clinical",
      diagnostic_confirmation == "Radiography without microscopic confirm"~ "Clinical",
      diagnostic_confirmation == "Unknown"~ "Unknown",
      TRUE ~ "Pathological"
  ),
  
  laterality = case_when(
    laterality == "Left - origin of primary"~ 1,
    laterality == "Right - origin of primary"~ 1,
    laterality == "Only one side - side unspecified"~ 1,
    laterality == "Not a paired site"~ 1,
    TRUE ~ 2
  ),
  
  cT = case_when(
    str_detect(derived_eod_2018_t_2018, "T1") ~ "T1",
    str_detect(derived_eod_2018_t_2018, "T2") ~ "T2",
    str_detect(derived_eod_2018_t_2018, "T3") ~ "T3",
    str_detect(derived_eod_2018_t_2018, "T4") ~ "T4",
    str_detect(derived_eod_2018_t_2018, "Tis") ~ "Tis",
    str_detect(derived_eod_2018_t_2018, "TX") ~ "TX"
    ),
  
  stage = case_when(
    str_detect(derived_eod_2018_stage_group_2018, "0") ~ "0",
    str_detect(derived_eod_2018_stage_group_2018, "1") ~ "I",
    str_detect(derived_eod_2018_stage_group_2018, "2") ~ "II",
    str_detect(derived_eod_2018_stage_group_2018, "3") ~ "III",
    str_detect(derived_eod_2018_stage_group_2018, "4") ~ "IV",
    TRUE ~ "Unknown"
    
  )) |> 
  select(sex, marital_status_at_diagnosis, 
         diagnostic_confirmation, grade_clinical_2018,
         laterality, combined_summary_stage_2004,
         derived_eod_2018_t_2018, derived_eod_2018_n_2018,
         derived_eod_2018_m_2018, derived_eod_2018_stage_group_2018,
         # rx_summ_scope_reg_ln_sur_2003, 
         rx_summ_surg_oth_reg_dis_2003,
         rx_summ_surg_rad_seq, rx_summ_systemic_sur_seq_2007,
         tumor_size_summary_2016:regional_nodes_positive_1988,
         seer_cause_specific_death_classification, survival_months,
         cs_status,
         cT, stage) |> na.omit()


    
   
miss <- surv |> 
  select(where(~any(is.na(.)))) |> names()


surv |> 
  gtsummary::tbl_summary(
    missing = "ifany"
    # type = list(everything() ~ "categorical"),
  ) |> 
  gtsummary::add_n()

# ─────────────────────────────────────────────────────────────────────────────
# 1. First we fit a Cox model
# ─────────────────────────────────────────────────────────────────────────────
library(survival)
library(glmnet)
library(ggsurvfit)

# z <- db[complete.cases(db), ] # For choosing cases without missed items#
# 
# I = sample(size =  round(nrow(z)/7,0),x =  1:nrow(z), replace = F) # Sampling from original data to construct test and train sets#
# 
# Datatrain = z[I,] #Introducing train set#
# Datatest = z[-I,] #Introducing test set#


# x <- as.matrix(surv |> select(-survival_months, -cs_status, -all_of(miss)))

x <- model.matrix( ~ ., surv)
y <- Surv(time = surv$survival_months, event = surv$cs_status)
fit <- glmnet(x, y, family = "cox", alpha = 1)

plot(fit)
coef(fit, s = 0.05)

survfit2(Surv(survival_months, cs_status) ~ cT, data = surv) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )

surv |>
  gtsummary::tbl_uvregression(
    method = survival::coxph,
    y = Surv(survival_months, cs_status),
    exponentiate = TRUE
  ) |> 
  gtsummary::bold_p() |> 
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Univariable Cox Regression Models",
                 subtitle = "Outcome: death (any cause)")
# ─────────────────────────────────────────────────────────────────────────────
# 1. Now we fit a RSF
# ─────────────────────────────────────────────────────────────────────────────
