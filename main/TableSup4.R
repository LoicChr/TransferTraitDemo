###########################################
#                                         #
#           Generation of Tab S4          #
#                                         #
###########################################
library(BayesianTools)
library(Rcpp)
library(ade4)
library(extraDistr)
library(stringr)
source("lib/trait2demo.R")
source("lib/likelihood.R")

perms <- c("123", "132", "213", "231","312", "321", "1234")
TabS4 <- sapply(perms, function(perm){
  if (perm == "123"){
    source("main/prior.R")
    chains <- list.files(paste0("results/obs/"), full.names = T, pattern = "chain")
  }else if(perm == "1234"){
    source("main/four_trait_analysis/prior_n4.R")
    chains <- list.files(paste0("results/obs_n4/"), full.names = T, pattern = "chain")[5:7]
    
  }else{
    source(paste0("main/permutation_analysis/prior_",perm,".R"))
    chains <- list.files(paste0("results/obs_", perm, "/"), full.names = T, pattern = "chain")
    
  }
  
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
  tr <- tr[,as.numeric(str_split(perm,"")[[1]])]
  
  if (perm == '1234'){
    source("lib/likelihood_n4.R")
  }else{
    source("lib/likelihood.R")
  }
  load(chains[1])
  nstep <- nrow(getSample(out))/3
  list_samplers <- lapply(chains, function(chain){
    load(chain)
    return(out)
  })
  Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')
  McmcSamplerList <- createMcmcSamplerList(list_samplers)
  
  # Posterior distribution
  Nchains <- length(list_samplers)
  burn <- nstep - 15000
  postDis <- getSample(McmcSamplerList, start = burn, thin= 250, parametersOnly = F)
  
  # Loglikelihood
  logLik_val <- median(postDis[,"Llikelihood"])
  
  # Convergence
  gelmrsf <- gelmanDiagnostics(McmcSamplerList, start = burn)$mpsrf
  
  Rcpp::sourceCpp('lib/LV_model_wrapped.cpp')
  # DIC 
  DIC.obs <- DIC(McmcSamplerList,start = burn)
  DIC.val <- round(DIC.obs$Dbar + DIC.obs$pV,1)
  
  # Pseudo R2
  H0 <- sum(apply(ixp, 1, function(x){
    dmultinom(x = x, prob = rep(1,ncol(ixp)), log = T)
  }))
  N <- sum(ixp)
  pseudoR2 <- median((1-exp(2/N*(H0-postDis[,"Llikelihood"])))/(1-exp(2/N*H0)))
  
  # Posterior parameters
  pars <- apply(postDis[,- which(colnames(postDis) %in% c("Lposterior","Llikelihood", "Lprior"))], 2, median)
  pars[grep("phi1", names(pars))] <- acos(pars[grep("phi1", names(pars))])
  if (perm == '1234'){
    pars[grep("phi2", names(pars))] <- acos(pars[grep("phi2", names(pars))])
  }
  
  coefs <- vector(mode = "numeric", length = 6*3)*NA
  names(coefs) <- apply(expand.grid(c("Tmin", "l", "c"), c("coef1", "coef2", "coef3", "coef4","a","b"), stringsAsFactors = F), 1, paste, collapse = "_")
  coefs <- coefs[-which(names(coefs) %in% c("c_b"))]  
  
  # Conversion in coefficients
  coefs[grep("coef1", names(coefs))] <- cos(pars[grep("phi1", names(pars))])
  coefs[grep("coef2", names(coefs))] <- cos(pars[grep("phi2", names(pars))])*sin(pars[grep("phi1", names(pars))])
  if (perm == '1234'){
    coefs[grep("coef3", names(coefs))] <- cos(pars[grep("phi3", names(pars))])*sin(pars[grep("phi2", names(pars))])*sin(pars[grep("phi1", names(pars))])
    coefs[grep("coef4", names(coefs))] <- sin(pars[grep("phi3", names(pars))])*sin(pars[grep("phi2", names(pars))])*sin(pars[grep("phi1", names(pars))])
    
  }else{
    coefs[grep("coef3", names(coefs))] <- sin(pars[grep("phi2", names(pars))])*sin(pars[grep("phi1", names(pars))])
  }
  coefs[grep("_a", names(coefs))] <- pars[grep("_a", names(pars))]
  coefs[grep("_b", names(coefs))] <- pars[grep("_b", names(pars))]
  
  # Reordering
  coefs_sp <- lapply(split(coefs, str_extract(names(coefs),"[:alpha:]+_")), function(coef_sp){
    coef_sp[1:3] <- coef_sp[match(c(1,2,3), as.numeric(str_split(perm,"")[[1]]))]
    coef_sp
  })
  coefs <- do.call(c,coefs_sp[3:1])
  names(coefs) <- unlist(sapply(coefs_sp[3:1], names))
  c(signif(coefs,3) , Nstep = nstep, Burn = burn, Nchains = Nchains, Gelman_msrf= signif(gelmrsf,3), logLik = signif(logLik_val,5), DIC = signif(DIC.val, 5), pseudoR2 = signif(pseudoR2, 2))
})
write.table(TabS4, file = "SourceData/TabS4.csv")
