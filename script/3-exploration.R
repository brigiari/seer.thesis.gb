library(tidyverse)
folder = "/sess-5"  

load(file = paste0(Sys.getenv("RAW_PATH"), folder, "/data.RData"))

source(here::here("script/4-recoding.R"))

# data[1:200,] |> 
#   select(-where(~all(is.na(.)))) %>% 
#   write.csv(paste0(Sys.getenv("RAW_PATH"), folder, "/test-2.csv"))
# 
# skimr::skim()
# 
# 
# table(data$primary_site_labeled)
# 
# data %>% names()
# 
# 
result <- restage(data)

result %>% 
  # select(-where(~all(is.na(.)))) %>%
  filter(!is.na(t_size_mm)) %>% 
  select(1:2,t_size_mm, derived_eod_2018_t_2018,  81:93) %>% 
  View()


result %>% 
  select()
  
