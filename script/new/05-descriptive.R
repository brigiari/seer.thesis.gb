# ─────────────────────────────────────────────────────────────────────────────
# 05-descriptive.R
# Descriptive tables, surgery-mortality analysis, missingness audit
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "04-analytic-dataset.R"))

cat("\n── Descriptive statistics ────────────────────────────────────────\n")

save_gt_table <- function(gt_obj, name, output_dir) {
  gt::gtsave(gt_obj, file.path(output_dir, paste0(name, ".html")))
  tryCatch(
    gt::gtsave(gt_obj, file.path(output_dir, paste0(name, ".png"))),
    error = function(err) {
      message("PNG export skipped for ", name, ": ", err$message)
    }
  )
}

# ── Table 1: Full cohort summary ────────────────────────────────────────────

tbl1 <- surv |>
  tbl_summary() |>
  bold_labels()

tbl1

# Save Table 1
tbl1_gt <- gtsummary::as_gt(tbl1)
save_gt_table(tbl1_gt, "table1", tbl_path)

# ── Surgery-mortality cross-table (MVV comment) ─────────────────────────────

cat("\n── Surgery vs Mortality ──────────────────────────────────────────\n")

surgery_mortality <- surv |>
  mutate(
    surgery_status   = if_else(surgery == "Yes", "Surgery", "No Surgery"),
    mortality_status = if_else(cs_status == 1, "Dead", "Alive")
  )

surgery_mortality_tbl <- surgery_mortality |>
  select(surgery_status, mortality_status) |>
  tbl_summary(by = surgery_status) |>
  add_overall() |>
  bold_labels()

surgery_mortality_tbl

# Save surgery-mortality table
surgery_tbl_gt <- gtsummary::as_gt(surgery_mortality_tbl)
save_gt_table(surgery_tbl_gt, "surgery-mortality", tbl_path)

# Surgery-mortality bar plot
surgery_mortality_plot <- surgery_mortality |>
  ggplot(aes(x = surgery_status, fill = factor(mortality_status,
                                               levels = c("Alive", "Dead")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  # Invert colors for Alive/Dead
  scale_fill_manual(values = {
    pal <- scales::hue_pal()(2)
    names(pal) <- c("Alive", "Dead")
    pal[c("Dead", "Alive")]
  }) +
  labs(
    title = "Mortality Rate by Surgery Status",
    x     = "Surgery Status",
    y     = "Proportion",
    fill  = "Status"
  ) +
  theme_minimal()

print(surgery_mortality_plot)

ggsave(file.path(fig_path, "surgery-mortality.png"),
       surgery_mortality_plot, width = 6, height = 4, dpi = 300)

# ── Pleura invasion missingness audit (MVV comment) ─────────────────────────

cat("\n── Pleura invasion missingness ───────────────────────────────────\n")

pleura_tbl <- surv |>
  count(visceral_and_parietal_pleural_invasion_recode_2010) |>
  mutate(pct = n / sum(n) * 100) |>
  arrange(desc(n))

print(pleura_tbl)

cat(sprintf(
  "  Unknown/NA in pleura invasion: %.1f%%\n",
  sum(pleura_tbl$pct[grepl("Unknown|Not applicable", 
                           pleura_tbl$visceral_and_parietal_pleural_invasion_recode_2010,
                           ignore.case = TRUE)])
))
