#############################################################################################
#                                                                                           #
#     Definition of the prior for the model with permuted trait axes in the order 1,3,2     #
#                                                                                           #
#############################################################################################

list_params <- apply(expand.grid(c("Tmin", "l", "c"), c("phi1", "phi2", "a", "b"), stringsAsFactors = F), 1, paste, collapse = "_")
list_params <- list_params[-which(list_params %in% c("c_b"))]  

# Parametrization of the prior for non-correlation parameters
Tmin_a_args <- c(mean = log(0.9), sd= 0.2)
l_a_args <- c(alpha = 9, beta = 3.5)
c_a_args <- c(alpha = 40, beta = 46)
Tmin_b_args <- c(mean = -2.37, sd= 0.35)
l_b_args <- c(mean = -6.5, sd= 0.35)

density = function(par){
  d1 = dunif(par[1], -1, 1, log =TRUE) #Tmin_phi1
  d2 = dunif(par[2], -1, 1, log =TRUE) #l_phi1
  d3 = dunif(par[3], -1, 1, log =TRUE) #c_phi1 
  
  d4 = dunif(par[4], pi, 3*pi, log =TRUE) #Tmin_phi2
  d5 = dunif(par[5], 0, 2*pi, log =TRUE) #l_phi2
  d6 = dunif(par[6], 0, 2*pi, log =TRUE) #c_phi2  
  
  d7 = dlnorm(par[7], Tmin_a_args["mean"], Tmin_a_args["sd"], log =TRUE) #Tmin_a
  d8 = dinvgamma(par[8], l_a_args["alpha"], beta = l_a_args["beta"], log =TRUE) #l_a 
  d9 = dinvgamma(par[9], c_a_args["alpha"], beta = c_a_args["beta"], log =TRUE) #c_a
  
  d10 = dnorm(par[10], Tmin_b_args["mean"], Tmin_b_args["sd"], log =TRUE) #Tmin_b
  d11 = dnorm(par[11], l_b_args["mean"], l_b_args["sd"], log =TRUE) #l_b 
  
  return(d1 + d2 + d3 + d4 +d5 +d6 +d7 +d8 +d9 +d10 +d11)
}
sampler = function(n=1){
  d1 = runif(n, -1, 1) #Tmin_phi1
  d2 = runif(n, -1, 1) #l_phi1
  d3 = runif(n, -1, 1) #c_phi1 
  
  d4 = runif(n, pi, 3*pi) #Tmin_phi2
  d5 = runif(n, 0, 2*pi) #l_phi2
  d6 = runif(n, 0, 2*pi) #c_phi2  
  
  d7 = rlnorm(n,Tmin_a_args["mean"], Tmin_a_args["sd"]) #Tmin_a
  d8 = rinvgamma(n, l_a_args["alpha"], beta = l_a_args["beta"]) #l_a 
  d9 = rinvgamma(n, c_a_args["alpha"], beta = c_a_args["beta"]) #c_a
  
  d10 = rnorm(n, Tmin_b_args["mean"], Tmin_b_args["sd"]) #Tmin_b
  d11 = rnorm(n, l_b_args["mean"], l_b_args["sd"]) #l_b 
  return(cbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11))
}

bounds[grep("phi1", list_params),1] <- -1
bounds[grep("phi1", list_params),2] <- 1
bounds[grep("phi2", list_params),1] <- 0
bounds[grep("phi2", list_params),2] <- 2*pi
bounds["Tmin_phi2",] <- bounds["Tmin_phi2",] + pi 

prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
