# ─────────────────────────────────────────────────────────────────────────────
# 00-setup.R
# Load libraries, set seed, configure parallelization, define paths
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(stringr)
library(gtsummary)
library(randomForestSRC)
library(survival)
library(ggplot2)
library(parallel)

seed <- 666
set.seed(seed)

options(rf.cores = detectCores(), mc.cores = detectCores())

data_folder <- "/data-raw/sess-10"
raw_path    <- Sys.getenv("RAW_PATH")

stopifnot(
  "RAW_PATH not set in .Renviron" = nchar(raw_path) > 0
)

fig_path   <- here::here("output", "figures")
tbl_path   <- here::here("output", "tables")

cat("Setup complete.\n")
cat("  RAW_PATH:", raw_path, "\n")
cat("  Cores:", detectCores(), "\n")
