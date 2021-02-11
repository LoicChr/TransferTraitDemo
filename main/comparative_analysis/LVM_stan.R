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

result = rstan::stan("lib/lvm.stan", 
                     data = list(E = ncol(X), SP = ncol(Y), L = 2 ,N = nrow(X), X = X, Y = Y), 
                     init = list(
                       list(W = matrix(rnorm(ncol(X)*ncol(Y), 0, 0.01), ncol(X), ncol(Y)),
                            LV = matrix(rnorm(nrow(X)*2, 0, 0.01), nrow(X), 2),
                            LF = matrix(rnorm(2*ncol(Y), 0, 0.01), 2, ncol(Y))
                            ),
                       list(W = matrix(rnorm(ncol(X)*ncol(Y), 0, 0.01), ncol(X), ncol(Y)),
                            LV = matrix(rnorm(nrow(X)*2, 0, 0.01), nrow(X), 2),
                            LF = matrix(rnorm(2*ncol(Y), 0, 0.01), 2, ncol(Y))
                       )
                     ),
                     chains = 2, cores = 2, iter = 100000, control = list(adapt_delta = 0.9))
saveRDS(result, "results/stan_lvm.RDS")