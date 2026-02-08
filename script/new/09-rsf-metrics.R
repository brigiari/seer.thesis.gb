# ─────────────────────────────────────────────────────────────────────────────
# 09-rsf-metrics.R
# Evaluate RSF: C-index, Brier score, CRPS, error rate plots
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "00-setup.R"))

train_path <- file.path(raw_path, "processed", "rsf-train-dev.RData")
if (file.exists(train_path)) {
  load(train_path)
  cat("Loaded cached RSF training results from:", train_path, "\n")
} else {
  source(here::here("script", "new", "08-rsf-train.R"))
}

obj <- rf_reduced

cat("\n── Model performance (OOB) ──────────────────────────────────────\n")

# OOB error rate (prediction error = 1 - C)
oob_err <- obj$err.rate[length(obj$err.rate)]
cat("  OOB error rate (PE):", round(oob_err, 4), "\n")

# C-index
c_index <- 1 - get.cindex(obj$yvar[, 1], obj$yvar[, 2], obj$predicted.oob)
cat("  C-index (OOB):     ", round(c_index, 4), "\n")

# ── Brier score ─────────────────────────────────────────────────────────────

cat("\n── Brier score ──────────────────────────────────────────────────\n")

bs_km  <- get.brier.survival(obj, cens.mode = "km")$brier.score
bs_rsf <- get.brier.survival(obj, cens.mode = "rfsrc")$brier.score

png(file.path(fig_path, "brier-score.png"), width = 8, height = 5,
    units = "in", res = 300)
plot(bs_km, type = "s", col = 2,
     xlab = "Time", ylab = "Brier Score",
     main = "Time-dependent Brier Score")
lines(bs_rsf, type = "s", col = 4)
legend("bottomright",
       legend = c("cens.model = km", "cens.model = rfsrc"),
       fill = c(2, 4))
dev.off()

# ── CRPS (integrated Brier score / time) ────────────────────────────────────

cat("\n── CRPS ─────────────────────────────────────────────────────────\n")

trapz <- randomForestSRC:::trapz
time  <- obj$time.interest

crps_km <- sapply(seq_along(time), function(j) {
  trapz(time[1:j], bs_km[1:j, 2] / diff(range(time[1:j])))
})
crps_rsf <- sapply(seq_along(time), function(j) {
  trapz(time[1:j], bs_rsf[1:j, 2] / diff(range(time[1:j])))
})

png(file.path(fig_path, "crps.png"), width = 8, height = 5,
    units = "in", res = 300)
plot(time, crps_km, ylab = "CRPS", type = "s", col = 2,
     main = "Continuous Rank Probability Score")
lines(time, crps_rsf, type = "s", col = 4)
legend("bottomright",
       legend = c("cens.model = km", "cens.model = rfsrc"),
       fill = c(2, 4))
dev.off()

# ── OOB survival plot ──────────────────────────────────────────────────────

png(file.path(fig_path, "oob-survival.png"), width = 8, height = 5,
    units = "in", res = 300)
plot.survival(obj)
dev.off()

# ── Error rate convergence ─────────────────────────────────────────────────

png(file.path(fig_path, "error-rate.png"), width = 8, height = 5,
    units = "in", res = 300)
plot(obj)
dev.off()

cat("\n  Figures saved to:", fig_path, "\n")

# ── Save workspace ─────────────────────────────────────────────────────────

save(list = ls(all.names = TRUE),
     file = paste0(raw_path, "/processed/rsf-dev-", Sys.Date(), ".RData"))
cat("  Workspace saved.\n")
