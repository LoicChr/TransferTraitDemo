############################################################################
#                                                                          #
#           Generation of figures 2-4 and associated source files          #
#                                                                          #
############################################################################
library(BayesianTools)
library(ade4)
library(purrr)
library(abind)
library(corrplot)
library(beeswarm)
library(tidyr)
library(RColorBrewer)

## Empirical functional traits
load("data/data.Rdata")
spxt_log <- spxt
spxt_log[,1:6] <-apply(spxt_log[,1:6], 2, log)
spxt_log2 <- apply(spxt_log, 2, function(x){
  x[is.na(x)] <- mean(x, na.rm = T)
  x
})
dudi.tr <- dudi.pca(spxt_log2, nf = 4, scannf = F)
tr <- apply(dudi.tr$li, 2, scale)

# load results
load("results/obs/pred.Rdata")

# Utility function to recompute demographic rates on a log scale
demoRates <- function(pars){
  source("lib/trait2demo.R")
  names(pars) <- colnames(postDis)
  pars[c("Tmin_phi1", "c_phi1", "l_phi1")] <- acos(pars[c("Tmin_phi1", "c_phi1", "l_phi1")])
  
  demo <- list()
  demo$l <- Trait2Demo(phi1 = pars["l_phi1"], phi2 = pars["l_phi2"],
                       a = pars["l_a"], b = pars["l_b"],
                       t1 = tr[,"Axis1"], t2 = tr[,"Axis2"], t3 = tr[,"Axis3"], exp = F)
  
  demo$c <- Trait2Demo(phi1 = pars["c_phi1"], phi2 = pars["c_phi2"],
                       a = pars["c_a"], b= -3.80,
                       t1 = tr[,"Axis1"], t2 = tr[,"Axis2"], t3 = tr[,"Axis3"], exp = F)
  
  demo$Tmin <- Trait2Demo(phi1 = pars["Tmin_phi1"], phi2 = pars["Tmin_phi2"],
                          a = pars["Tmin_a"], b = pars["Tmin_b"],
                          t1 = tr[,"Axis1"], t2 = tr[,"Axis2"], t3 = tr[,"Axis3"], exp = F)
  
  as.data.frame(demo)
}
# Correlations among demographic rates and between functional traits and demographic rates
cor_mats <- apply(postDis, 1, function(pars){
  #### Estimation of demographic rates
  demo <- demoRates(pars)
  
  cor_mat <- cor(spxt_log, demo, use = 'pairwise.complete.obs')
  demo_mat <- cor(demo)
  list(cor_trait = cor_mat, cor_demo = demo_mat)
})

# Figure 2A data
cor_trait_array <- do.call(abind, list(map(cor_mats, 1), along = 3))
cor_trait_tab <- do.call(rbind, lapply(map(cor_mats, 1), as.numeric))
colnames(cor_trait_tab) <- paste(rep(colnames(spxt),3), c("l", "c", "Tmin")[gl(3, ncol(spxt))], sep ="_")
write.csv(cor_trait_tab, file = "SourceData/Figure3A.csv", row.names = F)

cor_trait_array_med <- apply(cor_trait_array, c(1,2), quantile, probs = c(0.025,0.5, 0.975))

# Figure 3A data
cor_demo_array <- do.call(abind, list(map(cor_mats, 2), along = 3))
cor_demo_array_med <- apply(cor_demo_array, c(1,2), quantile, probs = c(0.025,0.5, 0.975), na.rm =T)
cor_demo_tab <- do.call(rbind, lapply(map(cor_mats, 2), function(x) as.numeric(x[upper.tri(x)])))
colnames(cor_demo_tab) <- c("l_c", "l_Tmin", "c_Tmin")
write.csv(cor_demo_tab, file = "SourceData/Figure3B.csv", row.names = F)


colPalette <- brewer.pal(11, "RdBu")[-c(1,2,10,11)]

# Figure 3 ----------------------------------------------------------------
cairo_ps("./figures/Figure3.eps", width = 9, height = 4.5)
layout(matrix(c(1,1,1,2,2), nrow = 1, byrow = T))
par(xpd = TRUE)
cor_med <- cor_trait_array_med[2,,]
colnames(cor_med) <- c("Minimum tolerated temperature", "Sensitivity to biomass", "Intraspecific competition rate")
row.names(cor_med) <- c("Reproductive height", "Vegetative height","Specific leaf area", "Leaf dry matter content", "Leaf carbon content", "Leaf nitrogen content", "$Leaf~delta^{13}~C", "$Leaf~delta^{15}~N")
corrplot(cor_med, method= "circle", cl.pos = "n",number.cex=0.9,mar = c(0,0,2,0), asp = 0.2, cl.cex = 0.9, cl.align.text = "l", tl.cex = 0.8, col = colPalette , addCoef.col = "black", tl.col = "black")
mtext("a", 2, adj=6, las=1, padj=-9, line = -2, cex = 0.8, font = 2)
par(xpd = TRUE)
cor_med <-cor_demo_array_med [2,,]
cor_med[upper.tri(cor_med, diag = T)] <- NA
colnames(cor_med) <- row.names(cor_med) <- c("Minimum tolerated temperature", "Sensitivity to biomass", "Intraspecific competition rate")
corrplot(cor_med, method= "circle", cl.pos = "n",na.label = " ",number.cex=0.9,mar = c(0,0,0,0), cl.cex =0.9, cl.align.text = "l", tl.cex = 0.8, col = colPalette , diag = T, addCoef.col = "black", tl.col = "black")
mtext("b", 2, adj=6, las=1, padj=-9, line = -2, cex = 0.8, font = 2)
dev.off()


# Figure 4 ----------------------------------------------------------------
pars.med <- apply(postDis, 2, median)
demo.med <- demoRates(pars.med)[,c("Tmin", "l", "c")]
demo.med$FG <- FG$FG
demo.med$species <- row.names(FG)
write.csv(demo.med, file = "SourceData/Figure4.csv", row.names = F)

pairwise.t.test(demo.med[,1] , FG$FG)
pairwise.t.test(demo.med[,2] , FG$FG)
pairwise.t.test(demo.med[,3] , FG$FG)
tests <- list(c('ab', 'a', 'a', 'b'), c('a', 'a', 'a', 'b'), c('a', 'b', 'b', 'ab'))

anova(lm(demo.med[,1] ~ FG$FG))
anova(lm(demo.med[,2] ~ FG$FG))
anova(lm(demo.med[,3] ~ FG$FG))


cairo_ps("./figures/Figure4.eps", width = 3.5, height = 10.5)
par(mfrow = c(3,1), mar = c(3.8,3,3,1), cex = 1., cex.main = 1., cex.axis = 0.8)
for (i in 1:3){
  boxplot(demo.med[,i] ~ FG$FG,las = 2, main = c("min. tol. Temperature","Sensitivity to biomass rate (log)", "Intraspecific competition rate (log)")[i], xlab = '',type ="n", ylim = c(min(demo.med[,i]),min(demo.med[,i])+diff(range(demo.med[,i]))*1.2), outline = FALSE, border = gray(0.25))
  beeswarm(demo.med[,i] ~ FG$FG, type = "square", pch = 21, add = T,col = "black", bg= c("#6699CC", "#CC6677", "#117733", "#DDCC77"), method = "center", cex = 0.9, xlab = '')
  text(c(1,2,3,4), y = min(demo.med[,i])+diff(range(demo.med[,i]))*1.05, labels = tests[[i]], pos = 3, cex = 1.)
  mtext(letters[i], 2, adj=6, las=1, padj=-10, line = -1, cex = 1.1, font = 2)
}
dev.off()


# Figure 2 ------------------------------------------------------
load("results/obs/pred.Rdata")
pred_uni <- as.data.frame(pred[seq(1,18, by = 2),])
temp_uni <- temp[seq(1,18, by = 2)]
colnames(pred_uni) <- row.names(spxt)
H0 <- apply(ixp, 1, function(x){
  dmultinom(x = x, prob = rep(1,ncol(ixp)), log = T)
})
# Data figure 2A
pred2 <- pred_uni
pred2$temp <- temp_uni
Fig2A_dat <- pivot_longer(data= as.data.frame(pred2), cols = contains("Sp"), names_to = "species")
colnames(Fig2A_dat) <- c("temp", "SpeciesCode", "relAb")
Fig2A_dat <- cbind(Fig2A_dat, FG[Fig2A_dat$SpeciesCode,])
write.csv(Fig2A_dat, file = "SourceData/Fig2A.csv", row.names = F)

# Data figure 2B
load("results/NM/Lls.Rdata")
R2.obs <-  (1-exp(2/rowSums(ixp)*(H0-logLik)))/(1-exp(2/rowSums(ixp)*(H0)))
NM.lls <- read.table("results/NM/Lls.txt", header = F)
Fig2B_dat <- data.frame(plot = names(temp), temp = temp, obs = R2.obs, NM_95 = apply(NM.lls, 2, quantile, probs = 0.95))
write.csv(Fig2B_dat, file = "SourceData/Fig2B.csv", row.names = F)


# Color palette
cols <-  c("#6699CC", "#CC6677", "#117733", "#DDCC77")

cairo_ps("figures/Figure2.eps", width = 7.5, height = 4)
par(mfrow = c(1,2), cex = 1., mar = c(4.5, 4.5, 2, 1), cex.main= 1, oma = c(1,1,1,1), cex.axis = 0.95)
plot(temp_uni, pred_uni[,1], type = "n", ylim = range(pred), main = "Modeled species relative abundance", xlab = "", ylab ="Relative abundance")
for (i in which(colSums(pred > 0.03) == 0)){ 
  x <- order(temp_uni)
  points(temp_uni[x], pred_uni[x,i], alpha = 0.8, type = "l", col = 'lightgray', lwd = 1)
}
for (i in which(colSums(pred > 0.03) > 0)){ 
  points(temp_uni[x], pred_uni[x,i], alpha = 0.8, type = "l",col = cols[FG$FG[i]], lwd = 1.2)
}
mtext("Mean annual temperature (°C)", side = 1, cex = 1, line = 3.5)
mtext("a", 2, adj=7, las=1, padj=-13, line =  -2, cex = 1.1, font = 2)
legend("topright", cex = 0.8,legend = levels(FG$FG), col = cols, lwd = 2.5)

x <-order(temp)
at.x <- c(1,2, 4,5, 7,8,10,11,13,14,16,17,19,20, 22,23,25,26)
bxp<-boxplot(lapply((1:ncol(NM.lls))[x], function(i) NM.lls[,i]), 
             outline = F, ylim = c(min(R2.obs),1), col = "white", border="white",
             xaxt ="n", ylab = expression(paste("Pseudo R"^"2")), at = at.x, main = "Model performance")
axis(side = 1, at = seq(1.5, 25.5,by = 3), labels = as.character(round(unique(temp)[order(unique(temp))],2)), las = 2)
mtext("Mean annual temperature (°C)", side = 1, cex = 1, line = 3.5)
points(at.x, R2.obs[x], bg = cols[2], pch = 23, cex = 1.2, lwd = 0.8)
points(at.x, apply(NM.lls, 2, quantile, probs = 0.95), bg = cols[1], pch = 24, cex = 1.2, lwd = 0.8)
abline(v = seq(3, 24,by = 3), lty = 2, col = "gray", lwd= 2)
abline(h = 0, lty = 4, lwd= 2)
mtext("b", 2, adj=7, las=1, padj=-13, line =  -2, cex = 1.1, font = 2)
dev.off()

