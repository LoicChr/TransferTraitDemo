##################################################################
#                                                                #
#     Definition of the prior for the model with four traits     #
#                                                                #
##################################################################

list_params <- apply(expand.grid(c("Tmin", "l", "c"), c("phi1", "phi2","phi3", "a", "b"), stringsAsFactors = F), 1, paste, collapse = "_")
list_params <- list_params[-which(list_params %in% c("c_b"))] 

# Parametrization of the prior for non-correlation parameters
Tmin_a_args <- c(mean = log(0.8), sd= 0.01)
l_a_args <- c(alpha = 40, beta = 12.3)
c_a_args <- c(alpha = 40, beta = 46)

Tmin_b_args <- c(mean = -2.37, sd= 0.01)
l_b_args <- c(mean = -6.1, sd= 0.2)

bounds <- data.frame(row.names = list_params, lower = rep(NA, length(list_params)), upper = NA)

bounds["Tmin_a", ] <- c(0.55,1.5)
bounds["l_a", ] <- c(0.2,0.74)
bounds["c_a", ] <- c(0.8,1.5)

bounds["Tmin_b", ] <- c(-3,-1.7)
bounds["l_b", ] <- c(-7.1,-5.5)

density = function(par){
  d1 = dtnorm(par[1], mean = -0.95, sd = 0.4, a = -1, b = 1, log =TRUE) #Tmin_phi1
  d2 = dtnorm(par[2], mean = 0.95, sd = 0.4, a = -1, b = 1, log =TRUE) #l_phi1
  d3 =dtnorm(par[3], mean = 0.58, sd = 0.3, a = -1, b = 1, log =TRUE) #c_phi1 
  
  d4 = dunif(par[4], -1,1, log =TRUE) #Tmin_phi2
  d5 = dunif(par[5], -1, 1, log =TRUE) #l_phi2
  d6 = dtnorm(par[6], mean = 0.98, sd = 0.4, a = -1, b = 1, log =TRUE) #c_phi2  
  
  d7 = dunif(par[7], 0, 2*pi, log =TRUE) #Tmin_phi3
  d8 = dunif(par[8], 0, 2*pi, log =TRUE) #l_phi3
  d9 = dunif(par[9], 0, 2*pi, log =TRUE) #c_phi3
  
  d10 = dlnorm(par[10], Tmin_a_args["mean"], Tmin_a_args["sd"], log =TRUE) #Tmin_a
  d11 = dinvgamma(par[11], l_a_args["alpha"], beta = l_a_args["beta"], log =TRUE) #l_a 
  d12 = dinvgamma(par[12], c_a_args["alpha"], beta = c_a_args["beta"], log =TRUE) #c_a
  
  d13 = dnorm(par[13], Tmin_b_args["mean"], Tmin_b_args["sd"], log =TRUE) #Tmin_b
  d14 = dnorm(par[14], l_b_args["mean"], l_b_args["sd"], log =TRUE) #l_b
  
  return(d1 + d2 + d3 + d4 +d5 +d6 +d7 +d8 +d9 +d10 +d11 + d12 + d13 + d14)
}
sampler = function(n=1){
  d1 = rtnorm(n, mean = -0.976, sd = 0.6, a = -1, b = 1) #Tmin_phi1
  d2 = rtnorm(n, mean = 0.95, sd = 0.6,a = -1, b = 1) #l_phi1
  d3 = rtnorm(n, mean = 0.58, sd = 0.5, a = -1, b = 1) #c_phi1 
  
  d4 = runif(n, -1, 1) #Tmin_phi2
  d5 = runif(n, -1, 1) #l_phi2
  d6 = rtnorm(n, mean = 0.95, sd = 0.6,a = -1, b = 1) #c_phi2 
  
  d7= runif(n, 0, 2*pi) #Tmin_phi3
  d8 = runif(n, 0, 2*pi) #l_phi3
  d9 = runif(n, 0, 2*pi) #c_phi3  
  
  d10 = rlnorm(n,Tmin_a_args["mean"], Tmin_a_args["sd"]) #Tmin_a
  d11 = rinvgamma(n, l_a_args["alpha"], beta = l_a_args["beta"]) #l_a 
  d12 = rinvgamma(n, c_a_args["alpha"], beta = c_a_args["beta"]) #c_a
  
  d13 = rnorm(n, Tmin_b_args["mean"], Tmin_b_args["sd"]) #Tmin_b
  d14 = rnorm(n, l_b_args["mean"], l_b_args["sd"]) #l_b  
  return(cbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14))
}

bounds[grep("phi1", list_params),1] <- -1
bounds[grep("phi1", list_params),2] <- 1
bounds[grep("phi2", list_params),1] <- -1
bounds[grep("phi2", list_params),2] <- 1
bounds[grep("phi3", list_params),1] <- 0
bounds[grep("phi3", list_params),2] <- 2*pi

prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
