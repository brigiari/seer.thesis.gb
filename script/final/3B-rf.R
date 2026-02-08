?get.brier.survival
plot.survival(obj)

## get brier score
km.brier = get.brier.survival(obj, cens.model = "km")$brier.score
rfsrc.brier = get.brier.survival(obj, cens.model = "rfsrc")$brier.score

plot(km.brier, type = "s", col = 2)
lines(rfsrc.brier, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))



# Now we run an analysis  under each of the four splitting rules.
# For each splitting rule we run 100 replications and record the mean and 
# standard deviation of the concordance error rate (as before ntree equals 1000):
splitrule <- c("logrank", "bs.gradien", "logrankscore")
nrep <- 100
err.rate <- matrix(0, 3, nrep)
names(err.rate) <- splitrule
ntree = 100L
# v.f <- as.formula("Survrsf(time,status) ~ .") #replaced with fml
for (j in 1:3) { 
  for (k in 1:nrep) { 
    err.rate[j,k] <-  randomForestSRC::rfsrc(fml, 
                                             data=surv,
                                             ntree=ntree, 
                                             splitrule=splitrule[j],
                                             seed=seed)$err.rate[ntree]
    } 
  } 

err.rate <- rbind( mean=apply(err.rate, 1, mean), 
                   std=apply(err.rate, 1, sd))
colnames(err.rate) <- splitrule
print(round(err.rate,3))


# We now consider the informativeness of each predictor under the best splitting rule
out <- randomForestSRC::rfsrc(fml, 
                              surv, 
                              ntree=ntree, 
                              splitrule = "logrank", 
                              seed=seed,
                              importance = "permute",
                              forest=T)
plot(out)

# The partial plots for the top six predictors are displayed here
plot.variable(out, partial=T, n.pred=6)

# We now consider the incremental effect of each predictor using a nested analysis:
# We sort predictors by their importance values and consider the nested sequence 
# of models starting with the top variable, followed by the model with the top 2 
# variables, then the model with the top three variables, and so on:
imp <- out$importance
pnames <- names(imp)
pnames.order <- pnames[rev(order(imp))]
n.pred <- length(pnames)
err <- rep(0, n.pred)

for (k in 1:n.pred){
  rsf.f <- "Surv(survival_months, cs_status)~"
  rsf.f <- as.formula(paste(rsf.f,
                            paste(pnames.order[1:k],collapse="+")))
  err[k] <- randomForestSRC::rfsrc(rsf.f, surv, ntree=ntree,
                                   splitrule="logrank")$err.rate[ntree]
}

imp.out <- as.data.frame(
  cbind(round(rev(sort(imp)),4),
        round(err,4),
        round(-diff(c(0.5,err)),4)),
  row.names=pnames.order)
colnames(imp.out)  <- c("Imp","Err","Drop Err")
print(imp.out)

# The first column is the importance value of a predictor in the full model. 
# The kth value in the second column is the error rate for the kth nested model, 
# while the kth value in the third column is the difference between the error 
# rate for the kth and (k − 1)th nested model, where the error rate for the null 
# model, k = 0, is 0.5. One can see that the top 3-4 predictors account for 
# much of the predictive power.


obj <- out
