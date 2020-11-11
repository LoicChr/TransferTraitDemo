
library(BayesianTools)
library(sjSDM)
load("TransferTraitDemo2/data/data.Rdata")
env <- scale(temp)

env <- as.matrix(env)
colnames(env)[1] <- "V1"

ixp <- as.matrix(ixp)


dev = torch$device("cuda:2")
ll = function(mu, Y, sigma) {
  mut = torch$tensor(mu, dtype=torch$float32, device=dev)$to(dev)
  YT = torch$tensor(Y, dtype = torch$float32, device=dev)$to(dev)
  sigmaT = torch$tensor(sigma, dtype = torch$float32, device=dev)$to(dev)
  noise = torch$randn(size=list(2000L, mut$shape[0], sigmaT$shape[1]), device=dev)
  E = torch$softmax(torch$tensordot(noise, sigmaT$t(), 1L)$add(mut), 2L)
  logprob = torch$distributions$Multinomial(probs = E)$log_prob(YT)
  maxlogprob = logprob$max(dim=0L)$values
  Eprob = logprob$sub(maxlogprob)$exp()$mean(dim = 0L)
  loss = Eprob$log()$neg()$sub(maxlogprob)
  return(loss$sum()$data$cpu()$numpy())
}


X = model.matrix(~poly(V1, degree = 2 ), data.frame(env))
Y = ixp


n = nrow(X)
e = ncol(X)
sp = ncol(Y)

latent = 4L

n_pars = e*sp + sp*latent


likelihood = function(pars) {
  W = matrix(pars[1:(e*sp)], e, sp)
  sigma = matrix(pars[-(1:(e*sp))], sp, latent)
  mu = X %*% W
  return(-ll(mu, Y, sigma))
}

sampler = function(n = 1) {
  return(t(sapply(1:n, function(i) c(rnorm(e*sp, 0, 0.2), rnorm(sp*latent, 0, 0.1)))))
}

density = function(par) {
  l1 = sum(dnorm(par[1:(e*sp)], 0, 0.2, log = TRUE)) + sum(dnorm(par[-(1:(e*sp))], 0, 0.1, log = TRUE))
  return(l1)
}

prior = createPrior(density = density, sampler = sampler)

bs = createBayesianSetup(likelihood = likelihood, 
                         
                         prior = prior,plotLower = rep(-1.0, n_pars), plotUpper = rep(1.0, n_pars), plotBest = rep(0.0, n_pars))

run = runMCMC(bs, settings = list(iterations = 100000L,burnin=80000L, sampler="DEzs", startValue = matrix(rnorm(n_pars*3, 0, 0.001), 3)))
saveRDS(run$chain, "mcmc_sjSDM.RDS")
chain = readRDS("mcmc_sjSDM.RDS")
### predicting #### 
predict = function(mu, sigma) {
  mut = torch$tensor(mu, dtype=torch$float32, device=dev)$to(dev)
  sigmaT = torch$tensor(sigma, dtype = torch$float32, device=dev)$to(dev)
  noise = torch$randn(size=list(2000L, mut$shape[0], sigmaT$shape[1]), device=dev)
  E = torch$softmax(torch$tensordot(noise, sigmaT$t(), 1L)$add(mut), 2L)
  E$mean(0L)
}

createPrediction = function(par) {
  W = matrix(pars[1:(e*sp)], e, sp)
  sigma = matrix(pars[-(1:(e*sp))], sp, latent)
  mu = X %*% W
  return(predict(mu, sigma)$cpu()$data$numpy())
}

samples = getSample(chain,numSamples = 2000L)

lls <- sapply(1:2000, function(i){
  pred = createPrediction( samples[i, 1:n_pars] )
  sum(sapply(1:nrow(ixp), function(k) dmultinom(x = ixp[k,], prob = unlist(pred[k,]+10e-5), log = T)))
})
saveRDS(lls, "lls_sjSDM.RDS")

mean(lls)
var(lls)
DIC <- var(-2 * lls)/2 -2 * mean(lls)

b <- readRDS(file = "results/lls_sjSDM.RDS")
c <- readRDS(file = "results/mcmc_sjSDM.RDS")
