###################################################################################################
#                                                                                                 #
#           Results and figures for the comparative analysis (Supplementary information)          #
#                                                                                                 #                                                                          #
###################################################################################################

predH0 <- readRDS("results//saturated.RDS")$pred
files <- c("results//jSDM_LVM.RDS","results//SDM.RDS","results//sjSDM_BT.RDS","results//sjSDM.RDS")

load("data/species_PA.Rdata")
library(pROC)
library(purrr)
library(rcartocolor)
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
roc.trf <- roc(as.numeric(bota > 0), as.numeric(pred[,!is.na(sp.sel)]))
roc.saturated <- roc(as.numeric(bota > 0), as.numeric(predH0[,!is.na(sp.sel)]))

auc.trf <- auc(roc.trf)
auc.sat <- auc(roc.saturated)
auc.sdms <- lapply(roc.sdms, auc)

### ROC curve graph
jpeg("figures/FigureS6.jpeg", width = 600, height = 600)
cols = carto_pal(6, 'Safe')[c(2,6,1,3:5)]
lwd.curve = 3.5
par(cex = 1.4, mar = c(5,5,2,2))
  plot(1-roc.trf$specificities, roc.trf$sensitivities, type = "l", col = cols[1], lwd = lwd.curve, xlab = "False positive rate", ylab = "True positive rate")
  points(1-roc.saturated$specificities, roc.saturated$sensitivities, type = "l", col = cols[2], lwd = lwd.curve)
  
  for (i in 1:4){
    points(1-roc.sdms[[i]]$specificities, roc.sdms[[i]]$sensitivities, type = "l", col = cols[i+2], lwd = lwd.curve)
  }
  abline(0,1)
  legend("bottomright", lwd = lwd.curve , cex = 0.9, col = cols, legend = c("Transfer function","Saturated model", "jSDM-LVM", "SDM", "sJSDM (MCMC)","sJSDM (MLE)"))
dev.off()
### Source data
FigS6_dat <- lapply(1:length(roc.sdms), function(i) data.frame(FPR = 1-roc.sdms[[i]]$specificities, TPR = roc.sdms[[i]]$sensitivities, Model = c("jSDM-LVM", "SDM", "sJSDM (MCMC)","sJSDM (MLE)")[i]))
FigS6_dat[length(roc.sdms)+1] <- data.frame(FPR = 1-roc.trf$specificities, TPR = roc.trf$sensitivities, Model = "TransferFunction")
FigS6_dat[length(roc.sdms)+2] <- data.frame(FPR = 1-roc.saturated$specificities, TPR = roc.saturated$sensitivities, Model = "SaturatedModel")

FigS6_dat <- do.call(rbind, FigS6_dat)
write.csv(FigS6_dat, file = "SourceData/FigS6.csv", row.names = F)


jpeg("figures/FigureS7.jpeg", width = 600, height = 600)
par(mar = c(5,4,2,4), oma = c(1,1,1,1), cex = 1.2, cex.lab = 1.3, cex.axis = 1.3)
lwd.curve = 5
plot(1:118, ylim = c(0,1), type = "n", xlab = "Number of species included", ylab = "Correlation with saturated model")
for (i in c(1:4)){
  points(out[[i]][[1]], type = "l", col = cols[i+2], lwd = lwd.curve )
}
points(cors_log, type = "l", col = cols[1], lwd = lwd.curve )
par(new = TRUE)
plot(colSums(ixp)[sp.ord], type = "l", col = "#888888", lty = 2, axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = lwd.curve, log ="y")
axis(side=4, at = c(1,5, 10,20,40,100), col = "#888888", col.ticks = "#888888", col.axis = "#888888")
mtext("Minimum number of sampled individuals", side=4, line=3, col = "#888888", cex = 1.4)
legend("topright", lwd = lwd.curve, cex = 1.2, col = cols[-2], legend = c("Transfer function","jSDM-LVM", "SDM", "sJSDM (MCMC)","sJSDM (MLE)"))
abline(v = 30, lty = 2)
dev.off()


FigS7_dat <- lapply(1:length(out), function(i) data.frame(NumberSpecies = 1:118, CorrelationSaturatedModel = out[[i]][[1]],Model = c("jSDM-LVM", "SDM", "sJSDM (MCMC)","sJSDM (MLE)")[i]))
FigS7_dat[length(out)+1] <- data.frame(NumberSpecies = 1:118, CorrelationSaturatedModel = cors_log, Model = "TransferFunction")
FigS7_dat <- data.frame(NumberSpecies = 1:118, 
           Correlation_TransferFunction=cors_log,
           Correlation_jSDMLVM = out[[1]][[1]],
           Correlation_SDM = out[[2]][[1]],
           Correlation_sJSDM_MCMC = out[[3]][[1]],
           Correlation_sJSDM_MLE = out[[4]][[1]],
           MinNumberIndividual = colSums(ixp)[sp.ord])
write.csv(FigS7_dat, file = "SourceData/FigS7.csv", row.names = F)

