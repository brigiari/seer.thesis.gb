library(tidyverse)
library(stringr)
library(gtsummary)
folder = "/data-raw/sess-10"  

load(file = paste0(Sys.getenv("RAW_PATH"), folder, "/data.RData"))


# 1. Initial filtering & basic recoding -----------------------------------------

db <- data |>
  # filter only primary tumor
  filter(primary_by_international_rules == "Yes") |>
  
  # Keep only lung & bronchus cases (double checked also with primary_site_labeled)
  filter(site_recode_icd_o_3_2023_revision == "Lung And Bronchus") |>  
  
  
#   Histologic codes in SEER: adenocarcinoma (AD) (8140, 8144, 8230, 8250, 8255, 8260, 8290, 8310, 8323, 
#                                                  8333, 8401, 8480, 8490, 8550, 8570, 8571, 8574), squamous cell carcinoma (SQCC) (8052, 8070–8075, 8083, 8084, 8123), large cell 
# carcinoma (LCC) (8012–8014, 8031–8033,8046,8082) and additional varieties of NSCLC (Other) (8022, 8200, 8240, 8430, 8560, 
#                                                                                             8562, 8980);
  filter(str_detect(histologic_type_icd_o_3, 
                    # adenocarcinoma-squamous cell carcinoma-large cell carcinoma-
                    "8140|8290|8310|8323|8333|8401|8480|8490|8550|8570|8571|8574| 
                    
                    8052|8070|8071|8072|8073|8074|8075|8083|8084|8123| 
                    
                    8012|8013|8014|8031|8032|8033|8046|8082|
                    
                    8022|8200|8240|8430|8560|8562|8980"
                    
                    )) |>  
  
  
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
    
    
    stage = if_else(
      stage %in% c("OCCULT", "UNK Stage"), NA_character_, stage
    ),  # convert "Unknown" to NA
    
    cT = if_else(cT=="T0", NA_character_, cT),
    
    
    regional_nodes_positive_1988 = case_when(
      regional_nodes_positive_1988 %in% c(1:90) ~ as.numeric(regional_nodes_positive_1988),
      regional_nodes_positive_1988 == 0 ~ 0,
      regional_nodes_positive_1988 %in% c(95, 97) ~ 1,
      regional_nodes_positive_1988 %in% c(98, 99) ~ NA_real_,
      TRUE ~ regional_nodes_positive_1988
    ),
    
    regional_nodes_examined_1988 = case_when(
      regional_nodes_examined_1988 %in% c(1:90) ~ as.numeric(regional_nodes_examined_1988),
      regional_nodes_examined_1988 == 0 ~ 0,
      regional_nodes_examined_1988 %in% c(95) ~ 0,
      regional_nodes_examined_1988 %in% c(96:99) ~ NA_real_,
      TRUE ~ regional_nodes_examined_1988
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
    
    radiation = if_else(rx_summ_surg_rad_seq == "No radiation and/or no surgery; unknown if surgery and/or radiation given",
                        "None", "Yes"),
    
    systemic_th = if_else(rx_summ_systemic_sur_seq_2007 == "No systemic therapy and/or surgical procedures",
                          "None", "Yes"),
    
    # re-factor any remaining character cols
    across(where(is.character), factor)
    
    
  ) |> 
  
  # ───────────────────────────────────────────────────────────────────────────
  # 5. Final cleanup & select modeling variables
  # ───────────────────────────────────────────────────────────────────────────
  
  drop_na(cT) |>       # drop T0 if present
  
  select(
    # core demographics & diagnosis
    age, sex, marital_status_at_diagnosis,
    diagnostic_confirmation, laterality, primary_site_labeled,

    
    # icd_o_3_hist_behav_malignant,
    primary_site_labeled,
    # treatment variables
    surgery, reg_lymph_surg, radiation, systemic_th,
    # cs variables
    cs_tumor_size_2004_2015:cs_mets_eval_2004_2015,
    # regional_nodes_positive,
    regional_nodes_examined_1988,
    regional_nodes_positive_1988,
    cs_site_specific_factor_1_2004_2017_varying_by_schema,
    cs_site_specific_factor_2_2004_2017_varying_by_schema,
  
    
    # outcome
    # seer_cause_specific_death_classification,
    survival_months, cs_status,
    
    # AJCC-derived TNM & stage
    cT, cN, M, stage
  )  
  # na.omit()              # drop any rows with remaining NAs


# table(db$rx_summ_systemic_sur_seq_2007)
# summary(surv$cs_mets_at_dx_2004_2015)


# ─────────────────────────────────────────────────────────────────────────────
# Explore data
# ─────────────────────────────────────────────────────────────────────────────

surv |> 
  gtsummary::tbl_summary(
    missing = "ifany"
    # type = list(everything() ~ "categorical"),
  ) |> 
  gtsummary::add_n()

# table(db$regional_nodes_examined_1988)

# ─────────────────────────────────────────────────────────────────────────────
# ORSF
# ─────────────────────────────────────────────────────────────────────────────

# install.packages("aorsf")
# library(aorsf)
# set.seed(666)
# fit_orfs <- orsf(data = surv,
#                  formula = survival_months + cs_status ~ .)
# 
# orsf_summarize_uni(fit_orfs, n_variables = 1)
# 
# 
# pd_by_gender <- orsf_pd_oob(fit_orsf, 
#                             pred_spec = list(sex = c("m", "f")),
#                             pred_horizon = 365 * 1:5)
# pd_by_gender %>% 
#   dplyr::select(pred_horizon, sex, mean) %>% 
#   tidyr::pivot_wider(names_from = sex, values_from = mean) %>% 
#   dplyr::mutate(ratio = m / f)


# data |> glimpse()


# # ─────────────────────────────────────────────────────────────────────────────
# # Fit a LASSO 
# # ─────────────────────────────────────────────────────────────────────────────
# library(survival)
# library(glmnet)
# library(ggsurvfit)
# # fit LASSO
# # 1. Identify the two columns to move up front
# key_cols <- c("survival_months", "cs_status")
# 
# # 2. Use setdiff() to get the rest of the column names
# other_cols <- setdiff(names(surv), key_cols)
# 
# x <- model.matrix( ~ ., surv[ , c(key_cols, other_cols)])
# y <- Surv(time = surv$survival_months, event = surv$cs_status)
# fit <- glmnet(x, y, family = "cox", alpha = 1)
# 
# 
# plot(fit)
# coef(fit, s = min(fit$lambda))
# 
# # final model
# fit <- glmnet(x, y, family = "cox", alpha = 1, 
#               lambda = min(fit$lambda))
# 
# coef(fit, s = min(fit$lambda))
# 
# 
# 
# # try with cv
# cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C")
# plot(cvfit)
# coef(cvfit, s = min(fit$lambda))
# 
# cvfit$lambda.min
# cvfit$lambda.1se
# 
# fit <- glmnet(x, y, family = "cox", alpha = 1, 
#               lambda = cvfit$lambda.min)
# 
# plot(fit)
# coef(fit, s = min(fit$lambda))













# ─────────────────────────────────────────────────────────────────────────────
# RSF
# ─────────────────────────────────────────────────────────────────────────────
library(randomForestSRC)

# Split into training and testing sets
seed <- 123  # for reproducibility

# Sample indices: 70% for training, 30% for testing
# NOTE: complete_data is not defined above—should be mydata or imputed dataset
ind <- sample(2, nrow(surv), replace = TRUE, prob = c(0.7, 0.3))
train <- surv[ind == 1, ]
test  <- surv[ind == 2, ]


# Training del modello -------------------------------------------------
tune <- randomForestSRC::tune(
  Surv(survival_months, cs_status) ~ .,
  data = train,
  ntreeTry = 200L,
  trace = TRUE,
  na.action = "na.impute",
  seed = seed
)

# Growing trees --------------------------------------------------
set.seed(seed)
rf_1 <- rfsrc(
  Surv(survival_months, cs_status) ~ .,
  data = train, 
  ntree = 100L,
  na.action = "na.impute",
  mtry = tune$optimal[["mtry"]],
  nodesize = tune$optimal[["nodesize"]],
  nodedepth = tune$rf$nodedepth,
  importance = "permute",
  seed = seed,
  do.trace = TRUE
  )

imp <- rf_1$importance


rf_complete <- rfsrc(
  Surv(survival_months, cs_status) ~ .,
  data = train, 
  ntree = 100L,
  mtry = tune$optimal[["mtry"]],
  nodesize = tune$optimal[["nodesize"]],
  nodedepth = tune$rf$nodedepth,
  seed = seed,
  do.trace = TRUE
)


library(randomForestSRC)
library(survival)
library(kernelshap)
library(shapviz)
# Continuous Rank Probability Scores
# Function that returns continuous rank probability scores
pred_fun <- function(model, data) {
  predict(model, data)$predicted
}
# Sample <=1000 rows from the training data. veteran is small enough to use all
xvars <- setdiff(colnames(surv), c("survival_months", "cs_status"))

X_explain <- surv[xvars]
sv <- kernelshap(rf_complete, X = X_explain, pred_fun = pred_fun) |> 
  shapviz()

sv |> sv_importance(kind = "bee")
sv |> sv_dependence(xvars)

# Performances ---------------------------------------------------------


# glimpse(surv)
# prediction error is measured by the C-error rate defined as 1-C, where C is 
# Harrell's (Harrell et al., 1982) concordance index.
v.obj <- rfsrc(Surv(survival_months, cs_status) ~ .,
               data = train, ntree=500, nodesize = 20, block.size = 10)

## print results of trained forest
print(v.obj)



## plot results of trained forest - higher the error, less robust the model
plot(v.obj)



# Variable Importance (VIMP) and Partial Plot
# VIMP (variable importance) is a technique for estimating the importance of a 
# variable by comparing performance of the estimated model with and without the 
# variable in it.

# # approach 1
# oo <- subsample(v.obj, verbose = FALSE)
# # take a delete-d-jackknife procedure for example
# vimpCI <- extract.subsample(oo)$var.jk.sel.Z
# vimpCI
# 
# 
# 
# # Confidence Intervals for VIMP
# plot.subsample(oo)
# # take the variable "XXX" for example for partial plot
# plot.variable(airq.obj, xvar.names = "XXX", partial = TRUE)


# VIMP
# most common measure is Breiman-Cutler VIMP, also called "permutation importance".
# NB: rather than using cross-validation, which can be computationally expensive, 
# permutation importance makes use of OOB estimation

# Large positive VIMP indicates high predictive ability while zero or negative 
# values identify noise variables. Subsampling [16] can be used to estimate the 
# standard error and to approximate the confidence intervals for VIMP. Figure 3 
# displays delete-d jackknife 99% asymptotic normal confidence intervals for the 
# p=39 variables from the systolic heart failure RSF analysis. Prediction error 
# was calculated using the C-index.


jk.obj <- subsample(v.obj)
pdf("VIMPsur.pdf", width = 15, height = 20)
par(oma = c(0.5, 10, 0.5, 0.5))
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
dev.off()

# ----------------------------------------------------------------------------
## plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)
## plot survival curves for first 10 individuals



## using function "plot.survival"
plot.survival(v.obj, subset = 1:10)
## obtain Brier score using KM and RSF censoring distribution estimators
bs.km <- get.brier.survival(v.obj, cens.model = "km")$brier.score
bs.rsf <- get.brier.survival(v.obj, cens.model = "rfsrc")$brier.score
## plot the brier score: measure used to assess prediction peformance.
# Lower values for the Brier score indicate better prediction performance
plot(bs.km, type = "s", col = 2)
lines(bs.rsf, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))


# Using the Brier score we can calculate the continuous rank probability 
# score (CRPS), defined as the integrated Brier score divided by time.
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
plot(tune.nodesize(Surv(survival_months, cs_status) ~ age +sex +stage 
              +systemic_th + radiation, 
              data = survtrace = TRUE)$err)



## ------------------------------------------------------------
## Minimal depth variable selection
## survival analysis
## use larger node size which is better for minimal depth
## ------------------------------------------------------------

# default call corresponds to minimal depth selection
vs.pbc <- var.select(object = v.obj)
topvars <- vs.pbc$topvars
# the above is equivalent to
max.subtree(v.obj)$topvars
# different levels of conservativeness
var.select(object = v.obj, conservative = "low")
var.select(object = v.obj, conservative = "medium")
var.select(object = v.obj, conservative = "high")




