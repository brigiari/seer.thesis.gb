# ----------------------------------------------
#  SEER Lung Cancer TNM Re‑staging Script
#  Editions: 7th (AJCC 2010), 8th (AJCC 2017), 9th (IASLC draft 2024)
#  Compatible with the column names listed by the user
#  Author: <your‑name>   Date: <date>
# ----------------------------------------------

# ===== 1 Packages =========================================================
library(tidyverse)

# ===== 2 Utility helpers ==================================================
na_blank <- function(x) {
  if (is.null(x)) return(x)
  y <- trimws(x)
  y[y == ""] <- NA_character_
  y
}
uc <- function(x) ifelse(is.na(x), NA_character_, toupper(x))

# ===== 3 Core recoding functions ==========================================

## 3.1 Tumour size (mm)
size_mm <- function(code) {
  case_when(
    code %in% 1:988           ~ as.numeric(code),
    code == 989               ~ 989,
    code == 990               ~ 5,
    code == 991               ~ 9,
    code == 992               ~ 15,
    code == 993               ~ 25,
    code == 994               ~ 35,
    code == 995               ~ 45,
    code %in% 996:998         ~ NA_real_,
    TRUE                      ~ NA_real_
  )
}

## 3.2 T category
classify_T7 <- function(mm) {
  out <- case_when(
    is.na(mm)         ~ NA_character_,
    mm <= 20          ~ "T1A",
    mm > 20 & mm <= 30 ~ "T1B",
    mm > 30 & mm <= 50 ~ "T2A",
    mm > 50 & mm <= 70 ~ "T2B",
    mm > 70           ~ "T3",
    TRUE              ~ NA_character_
  )
  uc(out)
}
classify_T8 <- function(mm) {
  out <- case_when(
    is.na(mm)         ~ NA_character_,
    mm <= 10          ~ "T1A",
    mm > 10 & mm <= 20 ~ "T1B",
    mm > 20 & mm <= 30 ~ "T1C",
    mm > 30 & mm <= 40 ~ "T2A",
    mm > 40 & mm <= 50 ~ "T2B",
    mm > 50 & mm <= 70 ~ "T3",
    mm > 70           ~ "T4",
    TRUE              ~ NA_character_
  )
  uc(out)
}
classify_T9 <- classify_T8   # unchanged from 8th

## 3.3 N category
map_N <- function(code) {
  out <- case_when(
    is.na(code)           ~ NA_character_,
    code %in% c(0,100)    ~ "N0",
    code %in% c(200,300)  ~ "N1",
    code %in% c(400,500)  ~ "N2",
    code %in% c(800,900)  ~ "N3",
    TRUE                  ~ NA_character_
  )
  uc(out)
}
classify_N9 <- function(ncat) {
  # assume single-station if unspecified
  ifelse(ncat == "N2", "N2A", ncat)
}

## 3.4 M category
map_M8 <- function(code) {
  out <- case_when(
    is.na(code)  ~ NA_character_,
    code == 0    ~ "M0",
    code == 10   ~ "M1A",
    code == 20   ~ "M1B",
    code == 30   ~ "M1C",
    TRUE         ~ NA_character_
  )
  uc(out)
}

# For 9th edition, with no metastasis count available, preserve M1C and vectorize
map_M9 <- function(code) {
  base <- map_M8(code)
  # vectorized fallback: keep M1C, otherwise pass through
  ifelse(base == "M1C", "M1C", base)
}

## 3.5 Stage group definition
make_stage <- function(T, N, M, edition = 8) {
  T <- uc(T); N <- uc(N); M <- uc(M)
  if (any(c(T,N,M) %in% c("TX","NX","MX", NA))) return(NA_character_)
  # metastasis → IV
  if (str_starts(M, "M1")) {
    return(ifelse(M %in% c("M1A","M1B"), "IVA", "IVB"))
  }
  key <- paste(T, N)
  stage_map <- list(
    "7" = list(
      "T1A N0"="IA","T1B N0"="IB","T2A N0"="IB","T2B N0"="IIA",
      "T1A N1"="IIA","T1B N1"="IIA","T2A N1"="IIB","T2B N1"="IIB",
      "T3 N0"="IIB","T3 N1"="IIIA","T1A N2"="IIIA","T2A N2"="IIIA",
      "T3 N2"="IIIB","T4 N0"="IIIA","T4 N1"="IIIA","T4 N2"="IIIB",
      "T1A N3"="IIIB","T2A N3"="IIIB","T3 N3"="IIIB","T4 N3"="IIIB"
    ),
    "8" = list(
      # N0
      "T1A N0"="IA1","T1B N0"="IA2","T1C N0"="IA3",
      "T2A N0"="IB","T2B N0"="IIA","T3 N0"="IIB","T4 N0"="IIIA",
      # N1
      "T1A N1"="IIB","T1B N1"="IIB","T1C N1"="IIB",
      "T2A N1"="IIB","T2B N1"="IIB","T3 N1"="IIIA","T4 N1"="IIIA",
      # N2
      "T1A N2"="IIIA","T1B N2"="IIIA","T1C N2"="IIIA",
      "T2A N2"="IIIA","T2B N2"="IIIB","T3 N2"="IIIB","T4 N2"="IIIB",
      # N3
      "T1A N3"="IIIC","T1B N3"="IIIC","T1C N3"="IIIC",
      "T2A N3"="IIIB","T2B N3"="IIIB","T3 N3"="IIIC","T4 N3"="IIIC"
    ),
    "9" = list(
      # N0
      "T1A N0"="IA1","T1B N0"="IA2","T1C N0"="IA3",
      "T2A N0"="IB","T2B N0"="IIA","T3 N0"="IIB","T4 N0"="IIIA",
      # N1
      "T1A N1"="IIA","T1B N1"="IIA","T1C N1"="IIA",
      "T2A N1"="IIA","T2B N1"="IIB","T3 N1"="IIIA","T4 N1"="IIIA",
      # N2A
      "T1A N2A"="IIB","T1B N2A"="IIB","T1C N2A"="IIB",
      "T2A N2A"="IIIA","T2B N2A"="IIIA","T3 N2A"="IIIA","T4 N2A"="IIIB",
      # N2B
      "T1A N2B"="IIIB","T1B N2B"="III A","T1C N2B"="IIIA",
      "T2A N2B"="IIIB","T2B N2B"="IIIB","T3 N2B"="IIIB","T4 N2B"="IIIB",
      # N3
      "T1A N3"="IIIB","T1B N3"="IIIB","T1C N3"="IIIB",
      "T2A N3"="IIIB","T2B N3"="IIIB","T3 N3"="IIIC","T4 N3"="IIIC"
    )
  )
  out <- stage_map[[as.character(edition)]][[key]]
  if (is.null(out)) return(NA_character_)
  out
}

# ===== 4 Restaging pipeline ==============================================
make_stage_vec <- Vectorize(make_stage, c("T","N","M"), SIMPLIFY=TRUE)
restage <- function(df) {
  df %>%
    mutate(across(matches("^(derived_|x7th_|cs_|eod_|_stage_group)"), na_blank)) %>%
    mutate(t_size_mm = coalesce(size_mm(cs_tumor_size_2004_2015),
                                size_mm(eod_primary_tumor_2018))) %>%
    mutate(
      T7 = coalesce(uc(derived_ajcc_t_7th_ed_2010_2015),
                    uc(derived_seer_combined_t_2016_2017),
                    classify_T7(t_size_mm)),
      N7 = coalesce(uc(derived_ajcc_n_7th_ed_2010_2015),
                    uc(derived_seer_combined_n_2016_2017),
                    map_N(cs_lymph_nodes_2004_2015)),
      M7 = coalesce(uc(derived_ajcc_m_7th_ed_2010_2015),
                    uc(derived_seer_combined_m_2016_2017),
                    map_M8(cs_mets_at_dx_2004_2015)),
      Stage7 = coalesce(uc(derived_ajcc_stage_group_7th_ed_2010_2015),
                        uc(x7th_edition_stage_group_recode_2016_2017),
                        make_stage_vec(uc(T7), uc(N7), uc(M7), 7))
    ) %>%
    mutate(
      T8 = coalesce(uc(derived_eod_2018_t_2018), classify_T8(t_size_mm)),
      N8 = coalesce(uc(derived_eod_2018_n_2018),
                    map_N(eod_regional_nodes_2018),
                    map_N(cs_lymph_nodes_2004_2015)),
      M8 = coalesce(uc(derived_eod_2018_m_2018), map_M8(eod_mets_2018),
                    map_M8(cs_mets_at_dx_2004_2015)),
      Stage8 = coalesce(uc(derived_eod_2018_stage_group_2018),
                        make_stage_vec(uc(T8), uc(N8), uc(M8), 8))
    ) %>%
    mutate(
      T9 = coalesce(uc(derived_eod_2018_t_2018), classify_T9(t_size_mm)),
      N9 = coalesce(uc(derived_eod_2018_n_2018),
                    classify_N9(map_N(eod_regional_nodes_2018))),
      M9 = coalesce(uc(derived_eod_2018_m_2018), map_M9(eod_mets_2018)),
      Stage9 = make_stage_vec(uc(T9), uc(N9), uc(M9), 9)
    )
}

# ===== 5 Example usage ====================================================
# result <- restage(data)
# write_csv(result, "seer_lung_tnm_restaged.csv")
# -------------------------------------------------------------------------
