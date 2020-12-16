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

result = rstan::stan("lib/glm.stan", 
                     data = list(E = ncol(X), SP = ncol(Y),N = nrow(X), X = X, Y = Y), 
                     init = list(
                       list(W = matrix(rnorm(ncol(X)*ncol(Y), 0, 0.01), ncol(X), ncol(Y))),
                       list(W = matrix(rnorm(ncol(X)*ncol(Y), 0, 0.01), ncol(X), ncol(Y)))
                     ),
                     chains = 2, cores = 2, iter = 100000)
saveRDS(result, "results/stan_glm.RDS")