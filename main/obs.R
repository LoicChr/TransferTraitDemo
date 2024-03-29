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
source("main/Data_prep.R")

#Id number for saving the output
id <- 1

# Reorganisation of the trait dataset
tr <- tr[,1:3]

######### Definition of the prior and likelihood
source("main/prior.R")
source("lib/likelihood.R")

bayesianSetup <- createBayesianSetup(LLpar, prior, names = list_params, parallel = 'external')

Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')

# # settings for the sampler
settings <- list(iterations = 50000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
save(list =   c("out"), file = paste0("results/obs/obs_chain", id, ".Rdata"))
