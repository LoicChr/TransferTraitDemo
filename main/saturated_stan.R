library(rstan)
library(sjSDM)
set.seed(42)
load("data/data.Rdata")
env <- scale(temp)

env <- as.matrix(env)
colnames(env)[1] <- "V1"

ixp <- as.matrix(ixp)
X = model.matrix(~poly(V1, degree = 2 ), data.frame(env))
Y = ixp

result = rstan::stan("lib/saturated.stan", 
                     data = list(SP = ncol(Y),N = nrow(X), Y = Y),
                     chains = 2, cores = 2, iter = 100000)
saveRDS(result, "results/stan_saturated.RDS")