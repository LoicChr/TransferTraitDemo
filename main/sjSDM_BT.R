library(BayesianTools)
library(sjSDM)
set.seed(42)
load("data/data.Rdata")
env <- scale(temp)

env <- as.matrix(env)
colnames(env)[1] <- "V1"

ixp <- as.matrix(ixp)
X = model.matrix(~poly(V1, degree = 2 ), data.frame(env))
Y = ixp

sampling = 2000L
model = sjSDM(Y = ixp, env = linear(data= env, formula = ~poly(V1, degree = 2 )), biotic = bioticStruct(df = 3L), se = FALSE, iter = 600L, family=multinomial(), sampling = 4000L, step_size = 1L, device = 2)

start_pars = c(as.vector(t(coef(model)[[1]])), as.vector(model$sigma))

dev = torch$device("cuda:2")
YT = torch$tensor(Y, dtype = torch$float32, device=dev)$to(dev)

ll = function(mu, Y, sigma) {
  mut = torch$tensor(mu, dtype=torch$float32, device=dev)$to(dev)
  sigmaT = torch$tensor(sigma, dtype = torch$float32, device=dev)$to(dev)
  noise = torch$randn(size=list(sampling, mut$shape[0], sigmaT$shape[1]), device=dev)
  E = torch$softmax(torch$tensordot(noise, sigmaT$t(), 1L)$add(mut), 2L)
  logprob = torch$distributions$Multinomial(probs = E)$log_prob(YT)
  maxlogprob = logprob$max(dim=0L)$values
  Eprob = logprob$sub(maxlogprob)$exp()$mean(dim = 0L)
  loss = Eprob$log()$neg()$sub(maxlogprob)
  return(loss$sum()$data$cpu()$numpy())
}


n = nrow(X)
e = ncol(X)
sp = ncol(Y)

latent = 3L
n_pars = e*sp + sp*latent


likelihood = function(pars) {
  W = matrix(pars[1:(e*sp)], e, sp)
  sigma = matrix(pars[-(1:(e*sp))], sp, latent)
  mu = X %*% W
  return(-ll(mu, Y, sigma))
}


sampler = function(n = 1) {
  return(sapply(start_pars, function(m) rnorm(n, m, 0.01)))
}

density = function(par) {
  return(sum(dnorm(par, 0, 3.0, log = TRUE)))
}

prior = createPrior(density = density, sampler = sampler)

bs = createBayesianSetup(likelihood = likelihood, prior = prior,plotLower = rep(-1.0, n_pars), plotUpper = rep(1.0, n_pars), plotBest = rep(0.0, n_pars))
run = runMCMC(bs, settings = list(iterations = 10.1e6,burnin=9.9e6, sampler="DEzs", Z = sampler(10L)))
saveRDS(run, "results/mcmc_sjSDM2.RDS")
