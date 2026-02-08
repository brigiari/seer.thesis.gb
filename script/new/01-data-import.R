# ─────────────────────────────────────────────────────────────────────────────
# 01-data-import.R
# Load raw SEER data from .RData file
# ─────────────────────────────────────────────────────────────────────────────

source(here::here("script", "new", "00-setup.R"))

load(file = paste0(raw_path, data_folder, "/datav2.RData"))

cat("Raw data loaded:", format(nrow(data), big.mark = ","), "records\n")
cat("Variables:", ncol(data), "\n")
