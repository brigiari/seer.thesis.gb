# source(here::here("script/2-tesi-data-cleaning.R"))
library(randomForestSRC)

source(here::here("script/final/2-tesi-data-cleaning.R"))

seed <- 666  # for reproducibility
surv |>
  gtsummary::tbl_summary()

# Training del modello -------------------------------------------------
tune <- randomForestSRC::tune(
  Surv(survival_months_scaled, cs_status) ~ .,
  data = surv,
  ntreeTry = 200L,
  trace = TRUE,
  na.action = "na.impute",
  seed = seed
)

# tune
# # $optimal
# # nodesize     mtry
# # 25       18
# #
# # $rf
# # NULL
# # Growing trees --------------------------------------------------
# set.seed(seed)
rf_1 <- randomForestSRC::rfsrc(
  Surv(survival_months_scaled, cs_status) ~ .,
  data = surv,
  ntree = 1000L,
  na.action = "na.impute",
  mtry = tune$optimal[["mtry"]],
  nodesize = tune$optimal[["nodesize"]],
  nodedepth = tune$rf$nodedepth,
  # importance = "permute",
  seed = seed,
  save.memory = TRUE,
  do.trace = TRUE
)
# xvar.used <- rf_1$var.used
# xvar.keep  <- names(xvar.used)[xvar.used >= 1]


# plot(rf_1)


# save(list = ls(all.names = TRUE),
#      file = paste0(Sys.getenv("RAW_PATH"), "/processed/full-2/nsclc-", Sys.Date(), ".RData"))


# print(rf_1)


# load(paste0(Sys.getenv("RAW_PATH"), "/processed/nsclc1000-whole.RData"))

# try var selection
imp_mat  <- vimp(rf_1, importance='permute')$importance         
imp_vec  <- setNames(as.vector(imp_mat),
                     names(imp_mat))       

imp_vec  <- imp_vec[!is.na(imp_vec)]          
imp_vec  <- imp_vec[imp_vec > 0]              
imp_rank <- sort(imp_vec, decreasing = TRUE)  
vars <- names(imp_rank[imp_rank > 0.001])
fml <- reformulate(vars, response = "Surv(survival_months_scaled, cs_status)") 


# set.seed(seed)
rf_1_var <- randomForestSRC::rfsrc(
  fml,
  data = surv,
  ntree = 1000L,
  na.action = "na.impute",
  mtry = tune$optimal[["mtry"]],
  nodesize = tune$optimal[["nodesize"]],
  nodedepth = tune$rf$nodedepth,
  # importance = "permute",
  seed = seed,
  save.memory = TRUE,
  do.trace = TRUE
)
print(rf_1_var)

save(list = ls(all.names = TRUE),
     file = paste0(Sys.getenv("RAW_PATH"), "/processed/full-2/nsclc-imp-", Sys.Date(), ".RData"))


obj <- rf_1_var

# Plot of survival estimates (OOB values)
plot.survival(obj)


# get metrics
# Prediction error is between 0 and 1, and measures how well the predictor 
# correctly ranks two random individuals in terms of survival. 
# Harrell’s C-index does not depend on choosing a fixed time for evaluation of 
# the model and specifically takes into account censoring of individuals

# The error rate is PE=1−C. Note that 0≤PE≤1 and that PE=0.5 corresponds to a 
# procedure doing no better than random guessing, whereas PE=0 indicates 
# perfect prediction.

# (OOB) Requested performance error is the error rate
obj$err.rate[length(obj$err.rate)]
# this is the built in function to get C-index
1-get.cindex(rf_1$yvar[,1], rf_1$yvar[,2], rf_1$predicted.oob)



# obj <- rf_1_var
# (OOB) CRPS is the integrated Brier score divided by time


bs.km <- get.brier.survival(obj, cens.mode = "km")$brier.score
bs.rsf <- get.brier.survival(obj, cens.mode = "rfsrc")$brier.score
# bs.test <- get.brier.survival(o.pred, cens.mode = "rfsrc")$brier.score

## plot the brier score
plot(bs.km, type = "s", col = 2)
lines(bs.rsf, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))



## here's how to calculate the CRPS for every time point
trapz <- randomForestSRC:::trapz
time <- obj$time.interest
crps.km <- sapply(seq_along(time), function(j) {
  trapz(time[1:j], bs.km[1:j, 2] / diff(range(time[1:j])))
})
crps.rsf <- sapply(seq_along(time), function(j) {
  trapz(time[1:j], bs.rsf[1:j, 2] / diff(range(time[1:j])))
})

## plot CRPS as function of time
plot(time, crps.km, ylab = "CRPS", type = "s", col = 2)
lines(time, crps.rsf, type ="s", col = 4)
legend("bottomright", legend=c("cens.model = rfsrc", "cens.model = km"), fill=c(4,2))



# # var impo
# jk.obj <- subsample(obj, B=2)
# o <- vimp(obj)
# 
# pdf("output/VIMPsur.pdf", width = 15, height = 20)
# par(oma = c(0.5, 10, 0.5, 0.5))
# par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
# plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
# dev.off()

# save(list = ls(all.names = TRUE),
#      file = paste0(Sys.getenv("RAW_PATH"), "/processed/nsclc800-final_workspace-full.RData"))
# load(paste0(Sys.getenv("RAW_PATH"), "/processed/nsclc800-final_workspace-full.RData"))

# Prediction ---------------------------------------------------------
# surv.pred <- predict(rf_1, surv, na.action = "na.impute")
# print(surv.pred)



# surv.pred.ran <- predict(rf_1, test, na.action = "na.random")
# print(surv.pred.ran)

# VIMP ------------------------------------------------------------
# importance <- predict(rf_1, 
#                       # get.tree=1:25,
#                       trace = TRUE,
#                       importance = TRUE)$importance

# 
# # this is for CI
# ## very small sample size so need largish subratio
# reg.smp.o <- subsample(rf_1 
#                        # B = 25, subratio = .5
# )
# ## summary of results
# print(reg.smp.o)
# plot.subsample(reg.smp.o)









library(randomForestSRC)
library(survival)
library(kernelshap)
library(shapviz)

# Continuous Rank Probability Scores




## ------------------------------------------------------------
## Minimal depth variable selection
## survival analysis
## use larger node size which is better for minimal depth
## ------------------------------------------------------------

# # default call corresponds to minimal depth selection
# vs.pbc <- var.select(object = obj)
# topvars <- vs.pbc$topvars
# # the above is equivalent to
# max.subtree(obj)$topvars
# # different levels of conservativeness
# var.select(object = obj, conservative = "low")
# var.select(object = obj, conservative = "medium")
# var.select(object = obj, conservative = "high")




## ------------------------------------------------------------
## Marginal plot
## ------------------------------------------------------------
plot.variable(rf_1, 
              xvar = "surgery", 
              time =  3*12, 
              surv.type = "surv", # Predicted survival (surv), where the predicted survival is for the time point specified using time (the default is the median follow up time).
              
              # cex.axis = 1.5, 
              # cex.lab = 1.5, 
              # cex.main = 1.5,
              # xlab = "Surgery",
              ylab = "Survival Probability"
              # main = "Surgery vs Survival Probability at 12 months"
              )  


## ------------------------------------------------------------
## Partial plot
## ------------------------------------------------------------
# plot.variable(rf_1, 
#               partial = TRUE,
#               xvar = "surgery", 
#               time.interest = c(1*12, 2*12, 3*12, 5*12, 10*12),
#               cex.axis = 1.5, 
#               cex.lab = 1.5, 
#               cex.main = 1.5,
#               xlab = "cT",
#               ylab = "Survival Probability",
#               main = "cT vs Survival Probability at 12 months")  

## ------------------------------------------------------------
## Now we see RF performance under complete-cases analysis
## ------------------------------------------------------------

set.seed(666)
surv.cc <- surv |> 
  na.omit()

trn.pt.cc <- sample(seq_len(nrow(surv.cc)), size = nrow(surv.cc)*0.7)
trn.cc <- surv.cc[trn.pt.cc, ]
tst.cc <- surv.cc[setdiff(seq_len(nrow(surv.cc)), trn.pt.cc), ]


# Split into training and testing sets
seed <- 666  # for reproducibility

# Training del modello -------------------------------------------------
tune.cc <- randomForestSRC::tune(
  Surv(survival_months_scaled, cs_status) ~ .,
  data = trn.cc,
  ntreeTry = 200L,
  trace = TRUE,
  na.action = "na.impute",
  seed = seed
)

# tune
# $optimal
# nodesize     mtry 
# 25       18 
# 
# $rf
# NULL
# Growing trees --------------------------------------------------
set.seed(seed)
rf.cc <- randomForestSRC::rfsrc(
  Surv(survival_months_scaled, cs_status) ~ .,
  data = trn.cc, 
  ntree = 1000L,
  na.action = "na.impute",
  mtry = tune.cc$optimal[["mtry"]],
  nodesize = tune.cc$optimal[["nodesize"]],
  nodedepth = tune.cc$rf$nodedepth,
  importance = "permute",
  seed = seed,
  save.memory = TRUE,
  do.trace = TRUE
)
plot(rf.cc)
plot.survival(rf.cc)




print(rf.cc)
o.pred.cc <- predict(rf.cc, tst.cc)
print(o.pred.cc)