# ----------------------------------------------
#  SEER Lung Cancer TNM Re‑staging Script
#  Editions: 7th (AJCC 2010), 8th (AJCC 2017), 9th (IASLC draft 2024)
#  Compatible with the column names listed by the user
#  Author: <your‑name>   Date: <date>
# ----------------------------------------------

# ===== 1 Packages =========================================================
library(tidyverse)   # dplyr, tidyr, readr …

# ===== 2 Utility helpers ==================================================

## Convert blanks / empty strings to proper NA
na_blank <- function(x) {
  if (is.null(x)) return(x)
  y <- trimws(x)
  y[y == ""] <- NA_character_
  y
}

## Always return upper‑case strings (keeps NA as NA)
uc <- function(x) ifelse(is.na(x), NA_character_, toupper(x))

# ===== 3 Core recoding functions ==========================================

## 3.1 Tumour size (mm) ----------------------------------------------------
size_mm <- function(code) {
  case_when(
    code %in% 1:988           ~ as.numeric(code),    # exact mm
    code == 989               ~ 989,                 # ≥989 mm
    code == 990               ~ 5,                   # microscopic focus
    code == 991               ~ 9,                   # "<1 cm" → 9 mm
    code == 992               ~ 15,                  # "<2 cm"
    code == 993               ~ 25,                  # "<3 cm"
    code == 994               ~ 35,                  # "<4 cm"
    code == 995               ~ 45,                  # "<5 cm"
    code %in% 996:998         ~ NA_real_,            # site‑specific / unknown
    TRUE                      ~ NA_real_
  )
}

## 3.2 T category ----------------------------------------------------------
classify_T7 <- function(mm) {
  out <- case_when(
    is.na(mm)         ~ NA_character_,
    mm <= 20          ~ "T1A",
    mm  > 20 & mm<=30 ~ "T1B",
    mm  > 30 & mm<=50 ~ "T2A",
    mm  > 50 & mm<=70 ~ "T2B",
    mm  > 70          ~ "T3",
    TRUE              ~ NA_character_
  )
  uc(out)
}

classify_T8 <- function(mm) {
  out <- case_when(
    is.na(mm)         ~ NA_character_,
    mm <= 10          ~ "T1A",
    mm  > 10 & mm<=20 ~ "T1B",
    mm  > 20 & mm<=30 ~ "T1C",
    mm  > 30 & mm<=40 ~ "T2A",
    mm  > 40 & mm<=50 ~ "T2B",
    mm  > 50 & mm<=70 ~ "T3",
    mm  > 70          ~ "T4",
    TRUE              ~ NA_character_
  )
  uc(out)
}

classify_T9 <- classify_T8   # identical size bands for 9th draft

## 3.3 N category ----------------------------------------------------------
map_N <- function(code) {
  out <- case_when(
    is.na(code)            ~ NA_character_,
    code %in% c(0, 100)    ~ "N0",
    code %in% c(200, 300)  ~ "N1",
    code %in% c(400, 500)  ~ "N2",
    code %in% c(800, 900)  ~ "N3",
    TRUE                   ~ NA_character_
  )
  uc(out)
}

classify_N9 <- function(ncat, multi_station_flag = NA) {
  ifelse(ncat != "N2", ncat,
         ifelse(is.na(multi_station_flag), "N2",
                ifelse(multi_station_flag == 1, "N2A", "N2B")))
}

## 3.4 M category ----------------------------------------------------------
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

map_M9 <- function(code, single_ext_flag = NA) {
  base <- map_M8(code)
  case_when(
    base != "M1B"             ~ base,
    is.na(single_ext_flag)     ~ "M1B",
    single_ext_flag == 1       ~ "M1C1",
    single_ext_flag  > 1       ~ "M1C2",
    TRUE                       ~ base
  )
}

## 3.5 Stage group (expanded TNM mapping for editions 7, 8, 9) -------------
make_stage <- function(T, N, M, edition = 8) {
  T <- uc(T); N <- uc(N); M <- uc(M)
  
  if (any(c(T, N, M) %in% c("TX", "NX", "MX", NA))) return (NA_character_)
  
  # --- M1 always stage IV ------------------------------------------------
  if (str_starts(M, "M1")) {
    return (ifelse(M %in% c("M1A", "M1B", "M1C1"), "IVA", "IVB"))
  }
  
  
  # --- Define comprehensive stage lookup by edition ----------------------
  key <- paste0(T, " ", N)
  stage_map <- list(
    "7" = list(
      "T1A N0" = "IA", "T1B N0" = "IB", "T2A N0" = "IB",
      "T2B N0" = "IIA", "T1A N1" = "IIA", "T1B N1" = "IIA",
      "T2A N1" = "IIB", "T2B N1" = "IIB", "T3 N0" = "IIB",
      "T3 N1" = "IIIA", "T1A N2" = "IIIA", "T2A N2" = "IIIA",
      "T3 N2" = "IIIB", "T4 N0" = "IIIA", "T4 N1" = "IIIA",
      "T4 N2" = "IIIB", "T1A N3" = "IIIB", "T2A N3" = "IIIB",
      "T3 N3" = "IIIB", "T4 N3" = "IIIB"
    ),
    "8" = list(
      "T1A N0" = "IA", "T1B N0" = "IA", "T1C N0" = "IA",
      "T2A N0" = "IB", "T2B N0" = "IIA", "T3 N0" = "IIB",
      "T4 N0" = "IIIA", "T1A N1" = "IIB", "T1B N1" = "IIB",
      "T2A N1" = "IIB", "T2B N1" = "IIB", "T3 N1" = "IIIA",
      "T4 N1" = "IIIA", "T1A N2" = "IIIA", "T3 N2" = "IIIB",
      "T4 N2" = "IIIB", "T1A N3" = "IIIC", "T4 N3" = "IIIC"
    ),
    "9" = list(
      "T1A N0" = "IA", "T1B N0" = "IA", "T1C N0" = "IA",
      "T2A N0" = "IB", "T2B N0" = "IIA", "T3 N0" = "IIB",
      "T4 N0" = "IIIA", "T1A N1" = "IIA", "T2A N1" = "IIA",
      "T2B N1" = "IIB", "T3 N1" = "IIIA", "T4 N1" = "IIIA",
      "T1A N2A" = "IIIA", "T1A N2B" = "IIIB", "T3 N2B" = "IIIB",
      "T4 N2A" = "IIIA", "T4 N2B" = "IIIB", "T1A N3" = "IIIC",
      "T4 N3" = "IIIC"
    )
  )
  
  edition_key <- as.character(edition)
  out <- stage_map[[edition_key]][[key]]
  if (is.null(out)) return(NA_character_)
  return(out)
  }




# ===== 4 Restaging pipeline ==============================================
## Vectorize `make_stage()` so it can be called row‑wise
make_stage_vec <- Vectorize(
  make_stage,
  vectorize.args = c("T", "N", "M"),
  SIMPLIFY = TRUE
)

## Updated `restage()` pipeline
restage <- function(df) {
  df %>%
    # blank out only truly missing helper fields
    mutate(across(matches("^(derived_|x7th_|cs_|eod_|_stage_group)"), na_blank)) %>%
    
    # compute tumor size once
    mutate(
      t_size_mm = coalesce(
        size_mm(cs_tumor_size_2004_2015),
        size_mm(eod_primary_tumor_2018)
      )
    ) %>%
    
    # compute 7th edition TNM and Stage
    mutate(
      T7 = coalesce(
        uc(derived_ajcc_t_7th_ed_2010_2015),
        uc(derived_seer_combined_t_2016_2017),
        classify_T7(t_size_mm)
      ),
      N7 = coalesce(
        uc(derived_ajcc_n_7th_ed_2010_2015),
        uc(derived_seer_combined_n_2016_2017),
        map_N(cs_lymph_nodes_2004_2015)
      ),
      M7 = coalesce(
        uc(derived_ajcc_m_7th_ed_2010_2015),
        uc(derived_seer_combined_m_2016_2017),
        map_M8(cs_mets_at_dx_2004_2015)
      ),
      Stage7 = coalesce(
        uc(derived_ajcc_stage_group_7th_ed_2010_2015),
        uc(x7th_edition_stage_group_recode_2016_2017),
        make_stage_vec(uc(T7), uc(N7), uc(M7), edition = 7)
      )
    ) %>%
    
    # compute 8th edition TNM and Stage
    mutate(
      T8 = coalesce(
        uc(derived_eod_2018_t_2018),
        classify_T8(t_size_mm)
      ),
      N8 = coalesce(
        uc(derived_eod_2018_n_2018),
        map_N(eod_regional_nodes_2018),
        map_N(cs_lymph_nodes_2004_2015)
      ),
      M8 = coalesce(
        uc(derived_eod_2018_m_2018),
        map_M8(eod_mets_2018),
        map_M8(cs_mets_at_dx_2004_2015)
      ),
      Stage8 = coalesce(
        uc(derived_eod_2018_stage_group_2018),
        make_stage_vec(uc(T8), uc(N8), uc(M8), edition = 8)
      )
    ) %>%
    
    # compute 9th edition TNM and Stage
    mutate(
      T9 = classify_T9(t_size_mm),
      N9 = classify_N9(map_N(eod_regional_nodes_2018)),
      M9 = map_M9(eod_mets_2018),
      Stage9 = make_stage_vec(uc(T9), uc(N9), uc(M9), edition = 9)
    )
}

# ===== 5 Example usage ====================================================
# result <- restage(data)
# View(result %>% select(t_size_mm, starts_with("T"), starts_with("N"), starts_with("M"), starts_with("Stage")))
# write_csv(result, "seer_lung_tnm_restaged.csv")
# -------------------------------------------------------------------------
