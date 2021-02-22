#########################################################################
#                                                                       #
#     Main script to fit the community model with functional traits     #
#     Permutation of trait axis in the order 2, 3, 1                    #
#                                                                       #
#########################################################################
setwd("./TransferTraitDemo/")

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
tr <- tr[,c(2,3,1)]
colnames(tr) <- colnames(dudi.tr$li)

######### Definition of the prior
source("main/permutation_analysis/prior_231.R")
source("lib/likelihood.R")

bayesianSetup <- createBayesianSetup(LLpar, prior, names = list_params, parallel = 'external')

Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')

# # settings for the sampler
settings <- list(iterations = 50000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
save(list = c("settings", "out", "bayesianSetup", "density", "sampler", "bounds","Tmin_a_args","l_a_args","c_a_args","Tmin_b_args" ,"l_b_args"), file = paste0("results/obs_231/obs_chain", id, ".Rdata"))

