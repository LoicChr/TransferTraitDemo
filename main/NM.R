################################################################################################
#                                                                                              #
#     Null model: main script to fit the community model with randomized functional traits     #
#                                                                                              #
################################################################################################
setwd("./TransferTraitDemo/")

library(BayesianTools)
library(Rcpp)
library(ade4)
library(extraDistr)
source("lib/trait2demo.R")
source("main/Data_prep.R")

#Id number for saving the output
id <- 1

# Randomization of trait values
tr <- tr[sample(1:nrow(tr)),1:3]

######### Definition of the prior
source("main/prior.R")
source("lib/likelihood.R")

bayesianSetup <- createBayesianSetup(LLpar, prior, names = list_params, parallel = 'external')

Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')

# # settings for the sampler
settings <- list(iterations = 50000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

pars.med <- apply(getSample(out, start = 35000, thin= 250, parametersOnly = F), 2, median)
ll.rand <- likelihood(pars.med, sum = F)
write.table(t(as.matrix(ll.rand)), file = "results/NM/Lls.txt", append = T, col.names = F, row.names = F)
