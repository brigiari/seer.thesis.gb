# ─────────────────────────────────────────────────────────────────────────────
# 07-rsf-tune.R
# Hyperparameter tuning for Random Survival Forest
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "06-sampling.R"))

cat("\n── RSF hyperparameter tuning (dev sample) ───────────────────────\n")

tune_path <- file.path(raw_path, "processed", "rsf-tune-dev.RData")
dir.create(dirname(tune_path), recursive = TRUE, showWarnings = FALSE)

if (file.exists(tune_path)) {
  load(tune_path)
  cat("Loaded cached tuning results from:", tune_path, "\n")
} else {
  tune_rsf <- tune(
    Surv(survival_months, cs_status) ~ .,
    data       = surv_dev,
    ntreeTry   = 200L,
    trace      = TRUE,
    na.action  = "na.impute",
    seed       = seed
  )

  save(list = c("tune_rsf", "surv_dev"), file = tune_path)
  cat("Saved tuning results to:", tune_path, "\n")
}

cat("\n  Optimal parameters:\n")
cat("    mtry:     ", tune_rsf$optimal[["mtry"]], "\n")
cat("    nodesize: ", tune_rsf$optimal[["nodesize"]], "\n")
