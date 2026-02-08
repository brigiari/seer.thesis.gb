tbl_chr <- surv |> 
  gtsummary::tbl_summary(
    missing = "ifany"
    # type = list(everything() ~ "categorical"),
  ) |> 
  gtsummary::add_n()

# Years at which we want status indicators
years <- c(1, 3, 5, 7)

# Loop over each time‐point and create a new column
for (y in years) {
  surv[[paste0("status_", y, "yr")]] <-
    with(surv,
      ifelse(
        survival_months <= y*12 & cs_status == 1,  # event happened on or before y
        1,
        ifelse(
          survival_months > y,                # still at risk past y (no event by y)
          0,
          NA_real_                 # censored before y → unknown
        )
      )
    )
}

# View results
tbl_survival <- surv %>% 
  select(contains("status"),
         -marital_status_at_diagnosis) %>% 
  gtsummary::tbl_summary()






set.seed(666)
trn.pt <- sample(1:nrow(surv), size = nrow(surv)*0.7)
trn <- surv[trn.pt, ]
tst <- surv[setdiff(1:nrow(surv), trn.pt), ]


trn %>% 
  mutate(group = "Train") %>% 
  rbind(tst %>% 
          mutate(group = "Test")) %>% 
  gtsummary::tbl_summary(
    missing = "ifany",
    by = group
    # type = list(everything() ~ "categorical"),
  ) |> 
  gtsummary::add_n() %>% 
  gtsummary::add_p() %>% 
  gtsummary::add_q("BH")
