library(tidyverse)
folder = "/data-raw/sess-5"  

load(file = paste0(Sys.getenv("RAW_PATH"), folder, "/data.RData"))


data$derived_seer_combined_t_2016_2017 <- str_replace_all(
  data$derived_seer_combined_t_2016_2017, "^c|^p", "T"
)
# source(here::here("script/4-recoding.R"))

result <- restage(data)

result %>% 
  # select(-where(~all(is.na(.)))) %>%
  filter(!is.na(t_size_mm)) %>% 
  select(1:2,
         
         # T7
         # derived_seer_combined_t_2016_2017,
         # derived_ajcc_t_7th_ed_2010_2015,
         # tumor_size_over_time_recode_1988,
         # eod_primary_tumor_2018,
         # cs_tumor_size_2004_2015,
         
         # N7
         # derived_ajcc_n_7th_ed_2010_2015,
         # derived_seer_combined_n_2016_2017,
         # cs_lymph_nodes_2004_2015,
         
         # M7
         derived_ajcc_m_7th_ed_2010_2015,
         derived_seer_combined_m_2016_2017,
         cs_mets_at_dx_2004_2015,
         # 
         # # STage 7
         # # derived_eod_2018_stage_group_2018,
         # # derived_ajcc_stage_group_7th_ed_2010_2015,
         # x7th_edition_stage_group_recode_2016_2017,
         # 
         # 
         # 
         # # T8
         # derived_eod_2018_t_2018,
         # # derived_ajcc_t_8th_ed_2018,
         # # derived_seer_combined_t_2016_2017,
         # # tumor_size_over_time_recode_1988,
         # # eod_primary_tumor_2018,
         # # cs_tumor_size_2004_2015,
         # 
         # # N8
         # derived_eod_2018_n_2018,
         # eod_regional_nodes_2018,
         # cs_lymph_nodes_2004_2015,
         # 
         # # M8
         # derived_eod_2018_m_2018,
         # eod_mets_2018,
         # cs_mets_at_dx_2004_2015,
         # 
         # # Stage 8
         # derived_eod_2018_stage_group_2018,
         # 
         # # T9
         # derived_eod_2018_t_2018,
         # 
         # # N9
         # derived_eod_2018_n_2018,
         # 
         # # M9
         # derived_eod_2018_m_2018,
         # eod_mets_2018,
         # 
         # # Stage 9
         
         
         
         
         
         
         
         
         
         81:93) %>% 
  View()


# T7
table(result$derived_seer_combined_t_2016_2017,result$T7)

# N7
table(result$derived_ajcc_n_7th_ed_2010_2015,result$N7)
 