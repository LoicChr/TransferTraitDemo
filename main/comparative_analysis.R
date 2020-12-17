
predH0 <- readRDS("results//saturated.RDS")$pred
files <- c("results//jSDM_LVM.RDS","results//SDM.RDS","results//sjSDM_BT.RDS","results//sjSDM.RDS")

load("data/species_PA.Rdata")
library(pROC)
library(purrr)

out = lapply(1:length(files), function(i){
  out.sdm = readRDS(files[i])
  #Pseudo R2
  H0 <- sum(apply(ixp, 1, function(x){
    dmultinom(x = x, prob = rep(1,ncol(ixp)), log = T)
  }))
  N <- sum(ixp)
  
  H1 <- sum(sapply(1:nrow(ixp), function(i){
    dmultinom(x = ixp[i,], prob = out.sdm$pred[i,], log = T)
  }))
  pseudoR2 <- (1-exp(2/N*(H0-H1)))/(1-exp(2/N*H0))
  
  cors_log <- numeric(length = 118)
  sp.ord <- order(colSums(ixp), decreasing = T)
  for (i in 1:118){
   cors_log[i] <- cor(log(as.numeric(predH0[,sp.ord][,1:i])+1e-5), log(as.numeric(out.sdm$pred[,sp.ord][,1:i])+1e-5))
  }
  
  #AUC
  roc.bota <- roc(as.numeric(bota > 0), as.numeric(out.sdm$pred[,!is.na(sp.sel)]))

  return(list(cors_log, pseudoR2, H1, roc.bota))
})

pseudoR2.sdms <- map(out, 2)
logLiks.sdms <- map(out, 3)
roc.sdms <- map(out, 4)

### Observed
load("results/obs/pred.Rdata")
cors_log <- numeric(length = 118)
sp.ord <- order(colSums(ixp), decreasing = T)
for (i in 1:118){
  cors_log[i] <- cor(log(as.numeric(predH0[,sp.ord][,1:i])+1e-5), log(as.numeric(pred[,sp.ord][,1:i])+1e-5))
}

### AUC scores
load("results/obs/pred.Rdata")
roc.obs <- roc(as.numeric(bota > 0), as.numeric(pred[,!is.na(sp.sel)]))
roc.saturated <- roc(as.numeric(bota > 0), as.numeric(predH0[,!is.na(sp.sel)]))

auc.obs <- auc(roc.obs)
auc.sat <- auc(roc.saturated)
auc.sdms <- lapply(roc.sdms, auc)

### ROC curve graph
png("FigS6.png", width = 600, height = 600)
lwd.curve = 3
par(cex = 1.4, mar = c(5,5,2,2))
  plot(1-roc.obs$specificities, roc.obs$sensitivities, type = "l", col = "forestgreen", lwd = lwd.curve, xlab = "False positive rate", ylab = "True positive rate")
  points(1-roc.saturated$specificities, roc.saturated$sensitivities, type = "l", col = "black", lwd = lwd.curve)
  
  for (i in 1:4){
    points(1-roc.sdms[[i]]$specificities, roc.sdms[[i]]$sensitivities, type = "l", col = c("blue","purple", "steelblue", "darkblue")[i], lwd = lwd.curve)
  }
  abline(0,1)
  legend("bottomright", lwd = 3, cex = 0.9, col = c("forestgreen","purple", "blue", "steelblue","darkblue", "black"), legend = c("Transfer function", "SDM", "jSDM-LVM", "sJSDM (MCMC)","sJSDM (MLE)","Saturated model"))
dev.off()

png("FigS7.png", width = 600, height = 600)
par(mar = c(5,4,2,4), oma = c(1,1,1,1), cex = 1.2, cex.lab = 1.3, cex.axis = 1.3)
plot(1:118, ylim = c(0,1), type = "n", xlab = "Number of species included", ylab = "Correlation with saturated model")
for (i in c(1:4)){
  points(out[[i]][[1]], type = "l", col = c("blue","purple", "steelblue", "darkblue")[i], lwd = 3)
}
points(cors_log, type = "l", col = "forestgreen", lwd = 3)
par(new = TRUE)
plot(colSums(ixp)[sp.ord], type = "l", col = "red", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3, log ="y")
axis(side=4, at = c(1,5, 10,20,40,100), col = "red", col.ticks = "red", col.axis = "red")
mtext("Minimum number of sampled individuals", side=4, line=3, col = "red", cex = 1.4)
legend("topright", lwd = 3, cex = 1.2, col = c("forestgreen","purple", "blue", "steelblue","darkblue"), legend = c("Transfer function", "SDM", "jSDM-LVM", "sJSDM (MCMC)","sJSDM (MLE)"))
abline(v = 30, lty = 2)
dev.off()




