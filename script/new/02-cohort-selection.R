# ─────────────────────────────────────────────────────────────────────────────
# 02-cohort-selection.R
# Apply all inclusion/exclusion filters with CONSORT-style logging
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "01-data-import.R"))

# ── ICD-O-3 histology codes ─────────────────────────────────────────────────

nsclc_codes <- sort(c(
  # Large-cell carcinomas
  8012, 8013, 8014, 8021, 8082,
  # Squamous cell carcinomas
  8070, 8071, 8072, 8073, 8078, 8083, 8084,
  # Adenocarcinomas
  8015, 8050, 8140, 8141, 8143, 8144, 8145, 8147, 8190, 8200, 8201, 8211,
  8230, 8250, 8251, 8252, 8253, 8254, 8255, 8256, 8257, 8260, 8265, 8290,
  8310, 8323, 8333, 8401, 8440, 8470, 8471, 8480, 8481, 8490, 8503, 8507,
  8550, 8551, 8570, 8571, 8572, 8574, 8576,
  # Other specific types
  8010, 8011, 8020, 8022, 8023, 8030, 8031, 8032, 8033, 8035, 8046, 8120,
  8240, 8241, 8243, 8244, 8245, 8246, 8249, 8430, 8525, 8560, 8562, 8575,
  8972, 8980, 8982
))

small_cell_codes <- c(8002, 8041, 8043, 8044, 8045)

nos_codes <- c(
  8046, 8003, 8004, 8022, 8030, 8031, 8032, 8033, 8035, 8120, 8200,
  8240, 8241, 8243, 8244, 8245, 8246, 8249, 8430, 8525, 8560, 8562, 8575
)

# ── CONSORT flow tracking ───────────────────────────────────────────────────

flow <- tibble(step = character(), n = integer(), excluded = integer())

log_step <- function(df, step_name) {
  n_now <- nrow(df)
  n_prev <- if (nrow(flow) == 0) n_now else tail(flow$n, 1)
  flow <<- bind_rows(flow, tibble(
    step     = step_name,
    n        = n_now,
    excluded = n_prev - n_now
  ))
  cat(sprintf("  %-55s N = %s (-%s)\n",
              step_name, format(n_now, big.mark = ","),
              format(n_prev - n_now, big.mark = ",")))
  df
}

# ── Type conversions (before filtering) ─────────────────────────────────────

cat("\n── Cohort selection ──────────────────────────────────────────────\n")

cohort <- data |>
  mutate(
    across(where(is.character), factor),
    time_from_diagnosis_to_treatment_in_days_recode =
      as.numeric(time_from_diagnosis_to_treatment_in_days_recode),
    survival_months = as.numeric(survival_months),
    total_number_of_in_situ_malignant_tumors_for_patient =
      as.numeric(total_number_of_in_situ_malignant_tumors_for_patient),
    tumor_size_over_time_recode_1988 =
      as.numeric(tumor_size_over_time_recode_1988),
    cs_status = if_else(
      seer_cause_specific_death_classification ==
        "Dead (attributable to this cancer dx)", 1, 0
    ),
    age = as.numeric(str_remove_all(age_recode_with_single_ages_and_90, " years"))
  ) |>
  log_step("Raw data after type conversion")

# ── Sequential filters ──────────────────────────────────────────────────────

cohort <- cohort |>
  filter(site_recode_icd_o_3_2023_revision == "Lung And Bronchus") |>
  log_step("Site recode: Lung and Bronchus")

cohort <- cohort |>
  filter(str_detect(primary_site_labeled, "34.0|34.1|34.2|34.3|34.8|34.9")) |>
  log_step("Primary site: C34.0-C34.3, C34.8, C34.9")

cohort <- cohort |>
  filter(histologic_type_icd_o_3 %in% nsclc_codes) |>
  log_step("Histology: NSCLC codes")

cohort <- cohort |>
  filter(!histologic_type_icd_o_3 %in% nos_codes) |>
  log_step("Exclude NOS histologies")

cohort <- cohort |>
  filter(primary_by_international_rules == "Yes") |>
  log_step("Primary tumors only")

cohort <- cohort |>
  filter(!type_of_reporting_source %in% c(
    "Death certificate only", "Autopsy only"
  )) |>
  log_step("Exclude death-cert-only / autopsy-only")

cohort <- cohort |>
  filter(age >= 18, age <= 80) |>
  log_step("Age 18-80")

cohort <- cohort |>
  filter(!derived_ajcc_t_7th_ed_2010_2015 %in% c(
    "T0", "T1NOS", "T2", "T2NOS", "TX"
  )) |>
  log_step("Exclude T0, T1NOS, T2/T2NOS, TX")

cohort <- cohort |>
  filter(!derived_ajcc_n_7th_ed_2010_2015 %in% c("NX")) |>
  log_step("Exclude NX")

cohort <- cohort |>
  filter(!derived_ajcc_m_7th_ed_2010_2015 %in% c("M1NOS", "M1")) |>
  log_step("Exclude M1/M1NOS")

cohort <- cohort |>
  filter(!derived_ajcc_stage_group_7th_ed_2010_2015 %in% c(
    "OCCULT", "UNK Stage"
  )) |>
  log_step("Exclude OCCULT / UNK Stage")

cohort <- cohort |>
  drop_na(
    derived_ajcc_t_7th_ed_2010_2015,
    derived_ajcc_n_7th_ed_2010_2015,
    derived_ajcc_m_7th_ed_2010_2015,
    derived_ajcc_stage_group_7th_ed_2010_2015
  ) |>
  log_step("Drop NA in AJCC staging")

cohort <- cohort |>
  drop_na(survival_months, cs_status) |>
  log_step("Drop NA in survival_months / cs_status")

cat("\n── CONSORT flow summary ──────────────────────────────────────────\n")
print(flow, n = Inf)
