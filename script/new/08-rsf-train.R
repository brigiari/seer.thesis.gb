# ─────────────────────────────────────────────────────────────────────────────
# 08-rsf-train.R
# Train RSF with tuned parameters, then variable selection via VIMP
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "07-rsf-tune.R"))

cat("\n── Training RSF (dev sample, 1000 trees) ─────────────────────────\n")

train_path <- file.path(raw_path, "processed", "rsf-train-dev.RData")
dir.create(dirname(train_path), recursive = TRUE, showWarnings = FALSE)

if (file.exists(train_path)) {
  load(train_path)
  cat("Loaded cached RSF training results from:", train_path, "\n")
} else {
  rf_full <- rfsrc(
    Surv(survival_months, cs_status) ~ .,
    data       = surv_dev,
    ntree      = 1000L,
    na.action  = "na.impute",
    mtry       = tune_rsf$optimal[["mtry"]],
    nodesize   = tune_rsf$optimal[["nodesize"]],
    seed       = seed,
    save.memory = TRUE,
    do.trace   = TRUE
  )

  cat("\n── Permutation VIMP ─────────────────────────────────────────────\n")

  imp_mat <- vimp(rf_full, importance = "permute")$importance
  imp_vec <- setNames(as.vector(imp_mat), names(imp_mat))
  imp_vec <- imp_vec[!is.na(imp_vec)]
  imp_vec <- imp_vec[imp_vec > 0]
  imp_rank <- sort(imp_vec, decreasing = TRUE)

  cat("  Variables with VIMP > 0.001:\n")
  print(imp_rank[imp_rank > 0.001])

  # ── Reduced model with selected variables ─────────────────────────────────
  vars_selected <- names(imp_rank[imp_rank > 0.001])
  fml_reduced   <- reformulate(vars_selected,
                               response = "Surv(survival_months, cs_status)")

  cat("\n── Training reduced RSF ─────────────────────────────────────────\n")
  cat("  Selected", length(vars_selected), "variables\n")

  rf_reduced <- rfsrc(
    fml_reduced,
    data       = surv_dev,
    ntree      = 1000L,
    na.action  = "na.impute",
    mtry       = tune_rsf$optimal[["mtry"]],
    nodesize   = tune_rsf$optimal[["nodesize"]],
    seed       = seed,
    save.memory = TRUE,
    do.trace   = TRUE
  )

  save(list = c("rf_full", "rf_reduced", "imp_rank", "vars_selected", "fml_reduced", "surv_dev"),
       file = train_path)
  cat("Saved RSF training results to:", train_path, "\n")
}
print(rf_reduced)
