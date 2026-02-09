# ─────────────────────────────────────────────────────────────────────────────
# 10-shap.R
# SHAP values via kernelshap + shapviz (dev model)
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "00-setup.R"))

library(kernelshap)
library(shapviz)

cat("\n── SHAP analysis (dev model) ───────────────────────────────────\n")

train_path <- file.path(raw_path, "processed", "rsf-train-dev.RData")
if (file.exists(train_path)) {
  load(train_path)
  cat("Loaded cached RSF training results from:", train_path, "\n")
} else {
  source(here::here("script", "new", "08-rsf-train.R"))
}

fig_path <- file.path(here::here("output"), "figures")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)

fit <- rf_reduced
xvars <- fit$xvar.names
sample_size <- min(1000L, nrow(surv_dev))
row_ids <- sample(seq_len(nrow(surv_dev)), size = sample_size)
X_explain <- as.data.frame(surv_dev[row_ids, xvars])

bg_size <- min(200L, nrow(surv_dev))
bg_ids <- sample(seq_len(nrow(surv_dev)), size = bg_size)
X_bg <- as.data.frame(surv_dev[bg_ids, xvars])

pred_fun <- function(model, data) {
  pred <- predict(model, data, na.action = "na.impute")$predicted
  if (is.matrix(pred) || is.array(pred)) {
    pred <- pred[, dim(pred)[2]]
  }
  as.numeric(pred)
}

sv <- kernelshap(fit, X = X_explain, bg_X = X_bg, pred_fun = pred_fun) |>
  shapviz()

shap_path <- file.path(raw_path, "processed", "rsf-shap-dev.RData")
save(list = c("sv", "row_ids", "xvars", "sample_size"), file = shap_path)
cat("Saved SHAP outputs to:", shap_path, "\n")

# ── SHAP plots ──────────────────────────────────────────────────────────────

beeswarm <- sv |>
  sv_importance(kind = "beeswarm") +
  theme_minimal()

ggsave(file.path(fig_path, "shap-beeswarm.png"), beeswarm,
       width = 9, height = 6, dpi = 300)

bar <- sv_importance(sv, show_numbers = TRUE) +
  theme_minimal()

ggsave(file.path(fig_path, "shap-bar.png"), bar,
       width = 8, height = 6, dpi = 300)

waterfall <- sv |>
  sv_waterfall(row_id = 1) +
  theme_minimal()

ggsave(file.path(fig_path, "shap-waterfall.png"), waterfall,
       width = 8, height = 6, dpi = 300)

force <- sv |>
  sv_force(row_id = 1)

ggsave(file.path(fig_path, "shap-force.png"), force,
       width = 9, height = 5, dpi = 300)

cat("  SHAP plots saved to:", fig_path, "\n")
