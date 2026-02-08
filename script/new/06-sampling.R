# ─────────────────────────────────────────────────────────────────────────────
# 06-sampling.R
# Create stratified development subsample (~10k) for faster iteration
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "04-analytic-dataset.R"))

cat("\n── Development subsample ─────────────────────────────────────────\n")

n_dev <- 10000L

if (nrow(surv) <= n_dev) {
  surv_dev <- surv
  cat("  Dataset smaller than", n_dev, "— using full dataset for development.\n")
} else {
  set.seed(seed)
  # stratified sample preserving event rate
  idx_event   <- which(surv$cs_status == 1)
  idx_noevent <- which(surv$cs_status == 0)

  event_rate <- length(idx_event) / nrow(surv)
  n_event    <- round(n_dev * event_rate)
  n_noevent  <- n_dev - n_event

  sampled_idx <- c(
    sample(idx_event,   min(n_event,   length(idx_event))),
    sample(idx_noevent, min(n_noevent, length(idx_noevent)))
  )

  surv_dev <- surv[sampled_idx, ]
  cat(sprintf("  Sampled %s records (event rate: %.1f%% vs full: %.1f%%)\n",
              format(nrow(surv_dev), big.mark = ","),
              mean(surv_dev$cs_status) * 100,
              mean(surv$cs_status) * 100))
}
