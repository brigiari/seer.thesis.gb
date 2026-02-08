library(dplyr)
library(stringr)
library(gtsummary)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Initial filtering & basic recoding
# ─────────────────────────────────────────────────────────────────────────────

db <- data |>
  # Keep only lung & bronchus primary diagnoses
  filter(site_recode_icd_o_3_2023_revision == "Lung And Bronchus") |>  
  
  # Remove any records missing T stage
  filter(!is.na(derived_ajcc_t_7th_ed_2010_2015)) |>    
  
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
  
  # ───────────────────────────────────────────────────────────────────────────
  # 2. Drop unwanted columns & rows with all NAs
  # ───────────────────────────────────────────────────────────────────────────
  
  select(
    -where(~ all(is.na(.))),               # drop columns entirely NA
    -primary_site,                         # drop redundant site codes
    -histologic_type_icd_o_3, 
    -icd_o_3_hist_behav,
    -survival_months_flag, 
    -patient_id,                           # drop identifier
    -age_recode_with_1_year_olds,          # drop alternate age fields
    -age_recode_with_1_year_olds_and_90,
    -age_recode_with_single_ages_and_85,
    -tumor_size_over_time_recode_1988,     # original size variable now numeric
    -contains("cod"),                      # any column with “cod” in its name
    -seer_registry_with_ca_and_ga_as_whole_states
  ) |> 
  drop_na(survival_months, cs_status) |>     # remove any remaining NAs in time/event
  filter(survival_months != 0)               # exclude zero-month follow-up


# ─────────────────────────────────────────────────────────────────────────────
# 3. Quick table of all variables (counts & missingness)
# ─────────────────────────────────────────────────────────────────────────────

db |>
  gtsummary::tbl_summary(
    missing = "ifany"   # show n (%) missing if any
  ) |>
  gtsummary::add_n()    # add total N column


# ─────────────────────────────────────────────────────────────────────────────
# 4. Further recoding for survival modeling
# ─────────────────────────────────────────────────────────────────────────────

surv <- db |> 
  
  mutate(
    # simplify marital status: married → “1”, all else → “0”
    marital_status_at_diagnosis = if_else(
      marital_status_at_diagnosis == "Married (including common law)", "1", "0"
    ),
    
    # collapse diagnostic confirmation into Clinical vs Pathological vs Unknown
    diagnostic_confirmation = case_when(
      diagnostic_confirmation %in% c(
        "Clinical diagnosis only",
        "Direct visualization without microscopic confirmation",
        "Clinical diagnosis only, direct visualization without microscopic confirmation",
        "Radiography without microscopic confirm"
      ) ~ "Clinical",
      diagnostic_confirmation == "Unknown" ~ "Unknown",
      TRUE ~ "Pathological"
    ),
    
    # unify laterality: any indication of side → “1”, else → “2”
    laterality = case_when(
      laterality %in% c(
        "Left - origin of primary",
        "Right - origin of primary",
        "Only one side - side unspecified",
        "Not a paired site"
      ) ~ "1",
      TRUE ~ "2"
    ),
    
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
    cT = if_else(cT=="T0", NA_character_, cT),
    # re-factor any remaining character cols
    across(where(is.character), factor)
  ) |> 
  
  # ───────────────────────────────────────────────────────────────────────────
  # 5. Final cleanup & select modeling variables
  # ───────────────────────────────────────────────────────────────────────────
  
  drop_na(cT) |>     # drop T0 if present
  
  select(
    # core demographics & diagnosis
    age, sex, marital_status_at_diagnosis,
    diagnostic_confirmation, laterality,
    
    # treatment variables
    rx_summ_surg_oth_reg_dis_2003,
    rx_summ_surg_rad_seq, 
    rx_summ_systemic_sur_seq_2007,
    
    # staging and evaluation metrics
    cs_tumor_size_2004_2015:cs_mets_eval_2004_2015,
    
    # outcome
    seer_cause_specific_death_classification,
    survival_months, cs_status,
    
    # AJCC-derived TNM & stage
    cT, cN, M, stage
  ) |> 
  na.omit()              # drop any rows with remaining NAs






# ─────────────────────────────────────────────────────────────────────────────
# Explore data
# ─────────────────────────────────────────────────────────────────────────────

surv |> 
  gtsummary::tbl_summary(
    missing = "ifany"
    # type = list(everything() ~ "categorical"),
  ) |> 
  gtsummary::add_n()

table(surv$cs_status)

# ─────────────────────────────────────────────────────────────────────────────
# ORSF
# ─────────────────────────────────────────────────────────────────────────────

# install.packages("aorsf")
library(aorsf)
set.seed(666)
fit_orfs <- orsf(data = surv,
                 formula = survival_months + cs_status ~ .)

orsf_summarize_uni(fit_orfs, n_variables = 1)


pd_by_gender <- orsf_pd_oob(fit_orsf, 
                            pred_spec = list(sex = c("m", "f")),
                            pred_horizon = 365 * 1:5)
pd_by_gender %>% 
  dplyr::select(pred_horizon, sex, mean) %>% 
  tidyr::pivot_wider(names_from = sex, values_from = mean) %>% 
  dplyr::mutate(ratio = m / f)




# ─────────────────────────────────────────────────────────────────────────────
# RSF
# ─────────────────────────────────────────────────────────────────────────────
library(randomForestSRC)
# glimpse(surv)
v.obj <- rfsrc(Surv(survival_months, cs_status) ~ ., data = surv, block.size = 1)
## plot tree number 3
plot(get.tree(v.obj, 3))
## print results of trained forest
print(v.obj)
## plot results of trained forest
plot(v.obj)
## plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)
## plot survival curves for first 10 individuals
## using function "plot.survival"
plot.survival(v.obj, subset = 1:10)
## obtain Brier score using KM and RSF censoring distribution estimators
bs.km <- get.brier.survival(v.obj, cens.model = "km")$brier.score
bs.rsf <- get.brier.survival(v.obj, cens.model = "rfsrc")$brier.score
## plot the brier score
plot(bs.km, type = "s", col = 2)
lines(bs.rsf, type ="s", col = 4)
legend("topright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
## plot CRPS (continuous rank probability score) as function of time
## here's how to calculate the CRPS for every time point
trapz <- randomForestSRC:::trapz
time <- v.obj$time.interes

crps.km <- sapply(1:length(time), function(j) {
  trapz(time[1:j], bs.km[1:j, 2] / diff(range(time[1:j])))
})
crps.rsf <- sapply(1:length(time), function(j) {
  trapz(time[1:j], bs.rsf[1:j, 2] / diff(range(time[1:j])))
})
plot(time, crps.km, ylab = "CRPS", type = "s", col = 2)
lines(time, crps.rsf, type ="s", col = 4)
legend("bottomright", legend=c("cens.model = km", "cens.model = rfsrc"), fill=c(2,4))

## fast nodesize optimization for veteran data
## optimal nodesize in survival is larger than other families
## see the function "tune" for more examples
tune.nodesize(Surv(time,status) ~ ., veteran)
## Primary biliary cirrhosis (PBC) of the liver
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)
print(pbc.obj)
## save.memory example for survival
## growing many deep trees creates memory issue without this option!
data(pbc, package = "randomForestSRC")
print(rfsrc(Surv(days, status) ~ ., pbc, splitrule = "random",
            ntree = 25000, nodesize = 1, save.memory = TRUE))
##------------------------------------------------------------
## trees can be plotted for any family
## see get.tree for details and more examples
##------------------------------------------------------------
## survival where factors have many levels
data(veteran, package = "randomForestSRC")
vd <- veteran
vd$celltype=factor(vd$celltype)
vd$diagtime=factor(vd$diagtime)
vd.obj <- rfsrc(Surv(time,status)~., vd, ntree = 100, nodesize = 5)
plot(get.tree(vd.obj, 3))

##------------------------------------------------------------
## example of imputation in survival analysis
##------------------------------------------------------------
data(pbc, package = "randomForestSRC")
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc, na.action = "na.impute")
## same as above but iterate the missing data algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc,
                  na.action = "na.impute", nimpute = 3)
## fast way to impute data (no inference is done)
## see impute for more details
pbc.imp <- impute(Surv(days, status) ~ ., pbc, splitrule = "random")
##------------------------------------------------------------
## compare RF-SRC to Cox regression
## Illustrates C-error and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------
if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE))
{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }
  ## data, formula specifications
  data(pbc, package = "randomForestSRC")
  pbc.na <- na.omit(pbc) ##remove NA's
  surv.f <- as.formula(Surv(days, status) ~ .)
  pec.f <- as.formula(Hist(days,status) ~ 1)
  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data = pbc.na, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f, pbc.na, ntree = 150)
  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbc.na, formula = pec.f,
                       splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)
  ## compute out-of-bag C-error for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, pbc.na)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b%%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(pbc.na$days[-train], pbc.na$status[-train], predict(cox.obj, pbc.na[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}