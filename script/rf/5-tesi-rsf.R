source(here::here("script/2-tesi-data-cleaning.R"))
library(randomForestSRC)

# Split into training and testing sets
seed <- 123  # for reproducibility

# Sample indices: 70% for training, 30% for testing
# NOTE: complete_data is not defined above—should be mydata or imputed dataset
ind <- sample(2, nrow(surv), replace = TRUE, prob = c(0.8, 0.2))
train <- surv[ind == 1, ]
test  <- surv[ind == 2, ]


# Training del modello -------------------------------------------------
tune <- randomForestSRC::tune(
  Surv(survival_months, cs_status) ~ .,
  data = train,
  ntreeTry = 200L,
  trace = TRUE,
  na.action = "na.impute",
  seed = seed
)

# Growing trees --------------------------------------------------
set.seed(seed)
rf_1 <- randomForestSRC::rfsrc(
  Surv(survival_months, cs_status) ~ .,
  data = train, 
  ntree = 100L,
  na.action = "na.impute",
  mtry = tune$optimal[["mtry"]],
  nodesize = tune$optimal[["nodesize"]],
  nodedepth = tune$rf$nodedepth,
  # importance = "permute",
  seed = seed,
  do.trace = TRUE
)

print(rf_1)

# Prediction ---------------------------------------------------------
surv.pred <- predict(rf_1, test, na.action = "na.impute")
print(surv.pred)



save(list = ls(all.names = TRUE), 
     file = here::here("processed/workspace.RData"))

# surv.pred.ran <- predict(rf_1, test, na.action = "na.random")
# print(surv.pred.ran)

# VIMP ------------------------------------------------------------
importance <- predict(rf_1, 
              # get.tree=1:25,
              trace = TRUE,
              importance = TRUE)$importance


save(list = ls(all.names = TRUE), 
     file = here::here("processed/workspace.RData"))






# library(randomForestSRC)
# library(survival)
# library(kernelshap)
# library(shapviz)
# 
# # Continuous Rank Probability Scores
# # Function that returns continuous rank probability scores
# pred_fun <- function(model, data) {
#   predict(model, data)$predicted
# }
# 
# # Sample <=1000 rows from the training data. veteran is small enough to use all
# X_explain <- veteran[xvars]
# sv <- kernelshap(fit, X = X_explain, pred_fun = pred_fun) |> 
#   shapviz()
# 
# sv |> sv_importance(kind = "bee")
# sv |> sv_dependence(xvars)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #
# # glimpse(surv)
# # prediction error is measured by the C-error rate defined as 1-C, where C is 
# # Harrell's (Harrell et al., 1982) concordance index.
# v.obj <- rfsrc(Surv(survival_months, cs_status) ~ .,
#                data = train, ntree=500, nodesize = 20, block.size = 10)
# 
# ## print results of trained forest
# print(v.obj)
# 
# 
# 
# ## plot results of trained forest - higher the error, less robust the model
# plot(v.obj)
# 
# 
# 
# # Variable Importance (VIMP) and Partial Plot
# # VIMP (variable importance) is a technique for estimating the importance of a 
# # variable by comparing performance of the estimated model with and without the 
# # variable in it.
# 
# # # approach 1
# # oo <- subsample(v.obj, verbose = FALSE)
# # # take a delete-d-jackknife procedure for example
# # vimpCI <- extract.subsample(oo)$var.jk.sel.Z
# # vimpCI
# # 
# # 
# # 
# # # Confidence Intervals for VIMP
# # plot.subsample(oo)
# # # take the variable "XXX" for example for partial plot
# # plot.variable(airq.obj, xvar.names = "XXX", partial = TRUE)
# 
# 
# # VIMP
# # most common measure is Breiman-Cutler VIMP, also called "permutation importance".
# # NB: rather than using cross-validation, which can be computationally expensive, 
# # permutation importance makes use of OOB estimation
# 
# # Large positive VIMP indicates high predictive ability while zero or negative 
# # values identify noise variables. Subsampling [16] can be used to estimate the 
# # standard error and to approximate the confidence intervals for VIMP. Figure 3 
# # displays delete-d jackknife 99% asymptotic normal confidence intervals for the 
# # p=39 variables from the systolic heart failure RSF analysis. Prediction error 
# # was calculated using the C-index.
# 
# 
# jk.obj <- subsample(v.obj)
# pdf("VIMPsur.pdf", width = 15, height = 20)
# par(oma = c(0.5, 10, 0.5, 0.5))
# par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
# plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
# dev.off()
# 
# # ----------------------------------------------------------------------------
# ## plot survival curves for first 10 individuals -- direct way
# matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
#         xlab = "Time", ylab = "Survival", type = "l", lty = 1)
# ## plot survival curves for first 10 individuals
# 
# 
# 
# ## using function "plot.survival"
# plot.survival(v.obj, subset = 1:10)
# ## obtain Brier score using KM and RSF censoring distribution estimators
# bs.km <- get.brier.survival(v.obj, cens.model = "km")$brier.score
# bs.rsf <- get.brier.survival(v.obj, cens.model = "rfsrc")$brier.score
# ## plot the brier score: measure used to assess prediction peformance.
# # Lower values for the Brier score indicate better prediction performance
# plot(bs.km, type = "s", col = 2)
# lines(bs.rsf, type ="s", col = 4)
# legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
# 
# 
# # Using the Brier score we can calculate the continuous rank probability 
# # score (CRPS), defined as the integrated Brier score divided by time.
# ## plot CRPS (continuous rank probability score) as function of time
# ## here's how to calculate the CRPS for every time point
# trapz <- randomForestSRC:::trapz
# time <- v.obj$time.interes
# crps.km <- sapply(1:length(time), function(j) {
#   trapz(time[1:j], bs.km[1:j, 2] / diff(range(time[1:j])))
# })
# crps.rsf <- sapply(1:length(time), function(j) {
#   trapz(time[1:j], bs.rsf[1:j, 2] / diff(range(time[1:j])))
# })
# plot(time, crps.km, ylab = "CRPS", type = "s", col = 2)
# 
# lines(time, crps.rsf, type ="s", col = 4)
# legend("bottomright", legend=c("cens.model = km", "cens.model = rfsrc"), fill=c(2,4))
# 
# ## fast nodesize optimization for veteran data
# ## optimal nodesize in survival is larger than other families
# ## see the function "tune" for more examples
# plot(tune.nodesize(Surv(survival_months, cs_status) ~ age +sex +stage 
#                    +systemic_th + radiation, 
#                    data = survtrace = TRUE)$err)
# 
# 
# 
# ## ------------------------------------------------------------
# ## Minimal depth variable selection
# ## survival analysis
# ## use larger node size which is better for minimal depth
# ## ------------------------------------------------------------
# 
# # default call corresponds to minimal depth selection
# vs.pbc <- var.select(object = v.obj)
# topvars <- vs.pbc$topvars
# # the above is equivalent to
# max.subtree(v.obj)$topvars
# # different levels of conservativeness
# var.select(object = v.obj, conservative = "low")
# var.select(object = v.obj, conservative = "medium")
# var.select(object = v.obj, conservative = "high")