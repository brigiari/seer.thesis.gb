# Feature importance = global view across all prediction
# While feature importance gives you a global overview, SHAP values dive deep, 
# helping you understand why, for example, a model decided a loan applicant is 
# risky this time, but not in another similar case.
#
# in machine learning, feature importance tells you which variables 
# (or “features”) in your dataset have the most influence on your model’s predictions.
#
# Now, feature importance can be calculated in different ways depending on the 
# type of model you’re working with. Here are a few common methods:
# 
# 1. Gini Importance (Used in Random Forests): This method ranks features based on
# how much they reduce uncertainty in decision trees.
# 2. Permutation Importance: Here’s a simple yet effective idea — randomly shuffle 
# each feature and see how much the model performance drops. If shuffling a feature 
# greatly impacts the performance, it’s an important one!
# 3. Coefficient-Based Feature Importance (For Linear Models): In linear models 
# like logistic regression, the magnitude of the coefficients tells you how much 
# influence a feature has.
#
# SHAP values
# 
# 1. Local and Global Interpretability: SHAP values allow you to explain not just the
# overall importance of a feature (like traditional methods) but also how that 
# feature impacts each individual prediction. That’s huge when you’re trying to 
# understand why a model made a specific decision.
# 2. Consistent and Unbiased: Unlike Gini importance, SHAP values don’t fall into the 
# same traps of bias from high-cardinality or correlated features. In fact, they 
# are designed to be consistent — if a feature becomes more important, its SHAP 
# value will always reflect that change.
# 2. Handling Interactions: If two features are working together to influence a 
# prediction, SHAP values capture this. For example: In a model predicting health 
# outcomes, perhaps both age and blood pressure are critical when considered together. 
# 3. SHAP values won’t just tell you these features are important — they’ll tell you 
# how they work in tandem to drive predictions.

library(randomForestSRC)
library(survival)
library(kernelshap)
library(shapviz)
library(ggplot2)
library(gtsummary)
library(tidyverse)

# this is the most comprehensive data (no splitting train/test)
# load(paste0(Sys.getenv("RAW_PATH"), "/processed/nsclc-shap-2025-08-042025-08-22", ".RData"))
xvars <- setdiff(colnames(surv), c("survival_months_scaled", "cs_status"))

fit <- rf_1_var

# Function that returns continuous rank probability scores
pred_fun <- function(model, data) {
  predict(model, data)$predicted
}

# Sample <=1000 rows from the training data. veteran is small enough to use all
set.seed(666)
sample <- sample(seq_len(nrow(surv)), size = 100L)
X_explain <- surv[sample, xvars]
sv <- kernelshap(fit, X = X_explain, pred_fun = pred_fun) |>
  shapviz()

save(list = ls(all.names = TRUE),
     file = paste0(Sys.getenv("RAW_PATH"), "/processed/nsclc-shap-2025-08-04", Sys.Date(), ".RData"))
# load(paste0(Sys.getenv("RAW_PATH"), "/processed/nsclc800-final_workspace-full.RData"))


# Compact SHAP analysis
# surv |> 
#   select(derived_ajcc_m_7th_ed_2010_2015,
#          surgery,
#          derived_ajcc_stage_group_7th_ed_2010_2015,
#          derived_ajcc_n_7th_ed_2010_2015,
#          chemotherapy_recode_yes_no_unk,
#          derived_ajcc_t_7th_ed_2010_2015,
#          sex,
#          grade_recode_thru_2017,
#          regional_nodes_positive,
#          reg_lymph_surg,
#          aya_site_recode_2020_revision,
#          age,
#          household_income,
#          primary_site_labeled,
#          visceral_and_parietal_pleural_invasion_recode_2010
#          ) |> 
#   names()

impA <- sv |> sv_importance(kind = "beeswarm")+
  theme_minimal() +
  scale_y_discrete(labels = function(x) str_to_title(gsub("_", " ", x)))

# ggsave(impA, filename = paste0(
#   "~/Library/CloudStorage/OneDrive-SharedLibraries-UnitofBiostatisticsEpidemiologyandPublicHealth/THO-STA IASLC - Documents/SEER/dissertation-gb/manuscript/manuscript"
#   , "/figures/impA.png"), width = 8, height = 6, dpi = 300)



impB <- sv_importance(sv, show_numbers = TRUE)+
  theme_minimal() +
  scale_y_discrete(labels = function(x) str_to_title(gsub("_", " ", x)))

# ggsave(impB, filename = paste0(
#   "~/Library/CloudStorage/OneDrive-SharedLibraries-UnitofBiostatisticsEpidemiologyandPublicHealth/THO-STA IASLC - Documents/SEER/dissertation-gb/manuscript/manuscript"
#   , "/figures/impB.png"), width = 8, height = 6, dpi = 300)

# sv |> sv_dependence(xvars)

# Decompose single predictions
single <- sv %>% sv_waterfall(row_id = 12) +
  theme(axis.text = element_text(size = 11)) +
  scale_y_discrete(labels = function(x) str_to_title(gsub("_", " ", x)))

# ggsave(single, filename = paste0(
#   "~/Library/CloudStorage/OneDrive-SharedLibraries-UnitofBiostatisticsEpidemiologyandPublicHealth/THO-STA IASLC - Documents/SEER/dissertation-gb/manuscript/manuscript"
#   , "/figures/singleA.png"), width = 8, height = 6, dpi = 300)

g <- sv %>% sv_force(row_id = 12) 
g[["data"]][["label"]] <- str_to_title(gsub("_", " ", g[["data"]][["label"]]))
# ggsave(g, filename = paste0(
#   "~/Library/CloudStorage/OneDrive-SharedLibraries-UnitofBiostatisticsEpidemiologyandPublicHealth/THO-STA IASLC - Documents/SEER/dissertation-gb/manuscript/manuscript"
#   , "/figures/singleB.png"), width = 8, height = 6, dpi = 300)



# Also multiple row_id can be passed: The SHAP values of the selected rows are 
# averaged and then plotted as aggregated SHAP values: The prediction profile for 
# surgery Yes:
wf <- sv %>% sv_waterfall(sv$X$surgery == "Yes") +
  theme(axis.text = element_text(size = 11)) +
  scale_y_discrete(labels = function(x) str_to_title(gsub("_", " ", x)))
  
# ggsave(wf, filename = paste0(
#     "~/Library/CloudStorage/OneDrive-SharedLibraries-UnitofBiostatisticsEpidemiologyandPublicHealth/THO-STA IASLC - Documents/SEER/dissertation-gb/manuscript/manuscript"
#     , "/figures/waterfall.png"), width = 8, height = 6, dpi = 300)

