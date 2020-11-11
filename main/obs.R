#########################################################################
#                                                                       #
#     Main script to fit the community model with functional traits     #
#                                                                       #
#########################################################################
setwd("./TransferTraitDemo/")

library(BayesianTools)
library(Rcpp)
library(ade4)
library(extraDistr)
source("lib/trait2demo.R")

#Id number for saving the output
id <- 1

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

######### Definition of the prior
source("main/prior.R")
source("lib/likelihood.R")

bayesianSetup <- createBayesianSetup(LLpar, prior, names = list_params, parallel = 'external')

Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')

# # settings for the sampler
settings <- list(iterations = 50000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
save(list = ls(), file = paste0(result_file, "/obs_chain", id, ".Rdata"))
-2300