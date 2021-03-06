###########################################
#                                         #
#           Result visualisation          #
#                                         #
###########################################
library(BayesianTools)
library(Rcpp)
library(ade4)
library(extraDistr)
source("lib/trait2demo.R")
source("lib/likelihood.R")
source("main/prior.R")

### Dataset
load("data/data.Rdata")
env <- scale(temp)

# Main trait axes
spxt_log <- spxt
spxt_log[,1:6] <-apply(spxt_log[,1:6], 2, log)
spxt_log <- apply(spxt_log, 2, function(x){
  x[is.na(x)] <- mean(x, na.rm = T)
  x
})
dudi.tr <- dudi.pca(spxt_log, nf = 3, scannf = F)
tr <- apply(dudi.tr$li, 2, scale)

chains <- list.files("results/obs/", full.names = T, pattern = "chain")
load(chains[1])
list_samplers <- lapply(chains, function(chain){
  load(chain)
  return(out)
})
Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')
McmcSamplerList <- createMcmcSamplerList(list_samplers)
plot(McmcSamplerList)

# Posterior distribution
postDis <- getSample(McmcSamplerList, start = 35000, thin= 250, parametersOnly = F)
write.csv(postDis, file = "SourceData/posterior.csv", row.names = F)


# Convergence
gelmanDiagnostics(McmcSamplerList, start = 35000)

#DIC
DIC.obs <- DIC(McmcSamplerList, start= 35000)
DIC.obs$Dbar + DIC.obs$pV

#Pseudo R2
H0 <- sum(apply(ixp, 1, function(x){
  dmultinom(x = x, prob = rep(1,ncol(ixp)), log = T)
}))
N <- sum(ixp)
pseudoR2 <- (1-exp(2/N*(H0-postDis[,"Llikelihood"])))/(1-exp(2/N*H0))
hist(pseudoR2)

# Computing the likelihood of the data for each individual plot.
source("lib/likelihood.R")
pars <- apply(postDis, 2, median)
pred <- likelihood(pars, pred.mat = T)
logLik <- likelihood(pars, sum = F)

save(list = c("pred", "postDis", "logLik"), file = "results/obs/pred.Rdata")
