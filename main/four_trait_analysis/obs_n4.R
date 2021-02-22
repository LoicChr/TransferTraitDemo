##############################################################################
#                                                                            #
#     Main script to fit the community model with four functional traits     #
#                                                                            #
##############################################################################
setwd(".TransferTraitDemo/")
  
rm(list = ls())
library(BayesianTools)
library(Rcpp)
library(ade4)
library(extraDistr)
source("lib/trait2demo.R")
source("main/Data_prep.R")

#Id number for saving the output
id <- 1

# Reorganisation of the trait dataset
tr <- tr[,1:4]

######### Definition of the prior
source("main/four_trait_analysis/prior_n4.R")
source("lib/likelihood_n4.R")

bayesianSetup <- createBayesianSetup(LLpar, prior, names = list_params, parallel = 'external')

Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')

# # settings for the sampler
settings <- list(iterations = 100000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
save(list = ls(), file = paste0("results/obs_n4/obs_chain", id, ".Rdata"))

