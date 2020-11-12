###########################################
#                                         #
#           Result visualisation          #
#                                         #
###########################################
library(BayesianTools)

chains <- list.files("results/obs/", full.names = T, pattern = "chain")
b=load(chains[1])
list_samplers <- lapply(chains, function(chain){
  load(chain)
  return(out)
})
Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')
McmcSamplerList <- createMcmcSamplerList(list_samplers)
plot(McmcSamplerList)

# Posterior distribution
postDis <- getSample(McmcSamplerList, start = 35000, thin= 50, parametersOnly = F)

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

# RMSE
pars <- apply(postDis, 2, median)
H0.pred <- matrix(1/ncol(ixp), nrow = nrow(ixp), ncol = ncol(ixp))

pred = likelihood(pars, pred.mat = T)
pred0 <- matrix(1/ncol(pred), ncol = ncol(pred), nrow = nrow(pred))

N = 2000
rmse <- numeric(N)
rmse0 <- numeric(N)

for (k in 1:N){
  rmse[k] <- sqrt(sum(sapply(1:nrow(pred), function(i){
    sim_i <- as.numeric(rmultinom(1,sum(ixp[i,]), pred[i,]))
    sum((sim_i - as.numeric(ixp[i,]))^2)
  }))/prod(dim(pred)))
}
for (k in 1:N){
  rmse0[k] <- sqrt(sum(sapply(1:nrow(pred), function(i){
    sim_i <- as.numeric(rmultinom(1,sum(ixp[i,]), pred0[i,]))
    sum((sim_i - as.numeric(ixp[i,]))^2)
  }))/prod(dim(pred)))
}
boxplot(rmse,rmse0, names = c("Transfer model", "H0")) 