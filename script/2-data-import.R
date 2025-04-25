folder = "/sess-5"  

data <- read.csv(paste0(Sys.getenv("RAW_PATH"), folder, "/export.csv"),
                 na.strings = c("Blank(s)", "NA", "")) |> 
  janitor::clean_names()

save(data, file = paste0(Sys.getenv("RAW_PATH"), folder, "/data.RData"))


renv::snapshot()
