folder = "/data-raw/sess-10"

data <- read.csv(paste0(Sys.getenv("RAW_PATH"), folder, "/export.csv"),
                 na.strings = c("Blank(s)", "NA", "", "N/A")) |> 
  janitor::clean_names()

save(data, file = paste0(Sys.getenv("RAW_PATH"), folder, "/datav2.RData"))


# renv::snapshot()
# usethis::edit_r_environ("project")

data |> 
  select(-where(~all(is.na(.)))) |> 
  names()
