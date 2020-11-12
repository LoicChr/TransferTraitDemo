###################################################
#                                                 #
#     Likelihood function for the calibration     #
#                                                 #
###################################################

likelihood<- function(pars, pred.mat = F){
  names(pars) <- list_params
  pars[c("Tmin_phi1", "c_phi1", "l_phi1")] <- acos(pars[c("Tmin_phi1", "c_phi1", "l_phi1")])
  
  # Create the matrice of theoretical functional traits
  demo <- list()
  demo$l <- Trait2Demo(phi1 = pars["l_phi1"], phi2 = pars["l_phi2"],
                           a = pars["l_a"], b = pars["l_b"],
                           t1 = tr[,"Axis1"], t2 = tr[,"Axis2"], t3 = tr[,"Axis3"], exp = T)
  
  demo$c <- Trait2Demo(phi1 = pars["c_phi1"], phi2 = pars["c_phi2"],
                           a = pars["c_a"], b= -3.80,
                           t1 = tr[,"Axis1"], t2 = tr[,"Axis2"], t3 = tr[,"Axis3"], exp = T)
  
  demo$Tmin <- Trait2Demo(phi1 = pars["Tmin_phi1"], phi2 = pars["Tmin_phi2"],
                               a = pars["Tmin_a"], b = pars["Tmin_b"],
                               t1 = tr[,"Axis1"], t2 = tr[,"Axis2"], t3 = tr[,"Axis3"], exp = F) ## Running of the model
  demo <- as.data.frame(demo)
  demo$g <- 10^-3.43
  N = nrow(demo)
  
  A <- demo[,"l"] %*% matrix(1, ncol = N)
  diag(A) <- diag(A) + demo[, "c"]
  t0 <- 0
  t1 <- 2e6
  dt <- 0.05 # dt doesn't seems to do anything as the integration is adaptive
  site <- as.numeric(unique(env))
  S <- length(site)
  logLik <- 0
  pred <- matrix(0, nrow = S, ncol = N)
  for (i in S:1){
    r <- as.numeric(demo[, "g"]*(site[i] - demo[,"Tmin"]))
    K = r/as.numeric(demo[, "l"]+ demo[,"c"])
    x <- r >0
    if (sum(x) > 1){
      init <- (runif(N)*K)[x]
      out <- integrateModel(init, r[x]*10^-1, A[x,x]*10^-1,t0, t1, dt)
      pred[i,x]<- unlist(out[nrow(out),-1]/sum(out[nrow(out),-1]))
      logLik <- logLik+ sum(sapply(which(env == site[i]), function(k)  dmultinom(x = ixp[k,], prob = unlist(pred[i,]+10e-5), log = T)))
      } else if (sum(x) == 1){                                          
      pred[i,x]<- 1
      logLik <- logLik+ sum(sapply(which(env == site[i]), function(k)  dmultinom(x = ixp[k,], prob = unlist(pred[i,]+10e-5), log = T)))
      }
    else{logLik <- -Inf; break}
  }
  if (pred.mat){
    pred <- pred[gl(9,2),]
    return(pred)
  }else{
    return(logLik)
  }
}

LLpar <- function(parsMat){
  parsMat <- matrix(parsMat, ncol = length(list_params), byrow = F)
  lls <- as.vector(apply(parsMat, 1, likelihood))
  return(lls)
}
