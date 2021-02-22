############################################################################
#                                                                          #
#     Hyperparameters to define the priors for the three traits model      #
#                                                                          #
############################################################################

list_params <- apply(expand.grid(c("Tmin", "l", "c"), c("phi1", "phi2", "a", "b"), stringsAsFactors = F), 1, paste, collapse = "_")
list_params <- list_params[-which(list_params %in% c("c_b"))]  

# Parametrization of the prior for non-correlation parameters
Tmin_a_args <- c(mean = log(0.9), sd= 0.2)
l_a_args <- c(alpha = 9, beta = 3.5)
c_a_args <- c(alpha = 40, beta = 46)
Tmin_b_args <- c(mean = -2.37, sd= 0.35)
l_b_args <- c(mean = -6.5, sd= 0.35)

bounds <- data.frame(row.names = list_params, lower = rep(NA, length(list_params)), upper = NA)
bounds["Tmin_a", ] <- c(0.55,1.5)
bounds["l_a", ] <- c(0.2,0.74)
bounds["c_a", ] <- c(0.8,1.5)
bounds["Tmin_b", ] <- c(-3,-1.7)
bounds["l_b", ] <- c(-7.1,-5.8)
