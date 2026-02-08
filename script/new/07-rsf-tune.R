# ─────────────────────────────────────────────────────────────────────────────
# 07-rsf-tune.R
# Hyperparameter tuning for Random Survival Forest
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "06-sampling.R"))

cat("\n── RSF hyperparameter tuning (dev sample) ───────────────────────\n")

tune_rsf <- tune(
  Surv(survival_months, cs_status) ~ .,
  data       = surv_dev,
  ntreeTry   = 200L,
  trace      = TRUE,
  na.action  = "na.impute",
  seed       = seed
)

cat("\n  Optimal parameters:\n")
cat("    mtry:     ", tune_rsf$optimal[["mtry"]], "\n")
cat("    nodesize: ", tune_rsf$optimal[["nodesize"]], "\n")
