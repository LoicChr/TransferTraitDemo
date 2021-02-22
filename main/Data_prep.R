###############################################
#                                             #
#     Help script: Preparation of the data    #
#                                             #
###############################################

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
