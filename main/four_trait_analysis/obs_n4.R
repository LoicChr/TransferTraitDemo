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


#Id number for saving the output
id <- as.numeric(commandArgs(trailingOnly = TRUE))

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
dudi.tr <- dudi.pca(spxt_log, nf = 4, scannf = F)
tr <- apply(dudi.tr$li, 2, scale)

######### Definition of the prior
source("main/four_trait_analysis/prior_n4.R")
source("lib/likelihood_n4.R")

bayesianSetup <- createBayesianSetup(LLpar, prior, names = list_params, parallel = 'external')

Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')

result_file <- "results_exp13/obs_n4_cons"
if (!file.exists(result_file)) dir.create(result_file, recursive = T)

# # settings for the sampler
settings <- list(iterations = 100000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
save(list = ls(), file = paste0("results/obs_n4/obs_chain", id, ".Rdata"))

