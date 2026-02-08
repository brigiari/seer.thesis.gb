x=-1
library(parallel)
detectCores()
options(rf.cores = x)
options(mc.cores = x)

options(rf.cores=detectCores(), mc.cores=detectCores())

imp_mat |> View()
names(imp_mat)
