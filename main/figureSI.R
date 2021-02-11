###################################
#                                 #
#           Figures S1-S2         #
#                                 #
###################################
library(ppcor)

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

species.temp <- apply(ixp, 2, function(x){
  weighted.mean(temp, x)
})
jpeg("./figures/FigureS1.jpeg", width = 4.5, height = 4.5, units = "in", res = 600)
par(mfrow = c(1,1), cex = 1., cex.main = 1., cex.axis = 0.8)
  plot(species.temp, dudi.tr$li[,1], xlab = "Species temperature niche (°C)", ylab = "PCA axis 1 score")
  abline(line(species.temp, dudi.tr$li[,1]))
dev.off()

reg <- lm(dudi.tr$li[,4] ~ relevel(FG[,2], ref = "Legumes"))
tests <- list(c('a', 'a', 'a', 'b'), c('a', 'b', 'b', 'b'), c('a', 'b', 'c', 'c'), c("a","ab","a","b"))

jpeg("./figures/FigureS2.jpg", width = 13.5, height = 4.5, units = "in", res = 600)
par(mfrow = c(1,4), mar = c(3.8,3,3,1), cex = 1., cex.main = 1., cex.axis = 0.8)
for (i in 1:4){
  boxplot(dudi.tr$li[,i] ~ FG[,2],las = 2, main = paste("PCA axis", i), xlab = '',type ="n", ylim = c(min(dudi.tr$li[,i]),min(dudi.tr$li[,i])+diff(range(dudi.tr$li[,i]))*1.2), outline = FALSE, border = gray(0.25))
  beeswarm(dudi.tr$li[,i] ~ FG[,2], type = "square", pch = 21, add = T,col = "black", bg= c("lightblue", "lightpink", "lightgreen", "peru"), method = "center", cex = 0.9, xlab = '')
  text(c(1,2,3,4), y = min(dudi.tr$li[,i])+diff(range(dudi.tr$li[,i]))*1.05, labels = tests[[i]], pos = 3, cex = 1.)
}
dev.off()


#SourceData
FigS1S2_dat <- data.frame(SpeciesCode = row.names(dudi.tr$li), dudi.tr$li[,1:4], FG, TempNiche = species.temp)
write.csv(FigS1S2_dat, file = "SourceData/FigS1S2.csv", row.names = F)

TabS1_dat <- data.frame(Trait = row.names(dudi.tr$co), dudi.tr$co[,1:4])
colnames(TabS1_dat)[2:5] <- paste0('Axis',1:4)
TabS1_dat<-rbind(TabS1_dat, c("Explained variance", dudi.tr$eig[1:4]))
write.csv(TabS1_dat, file = "SourceData/TabS1.csv", row.names = F)

### Figure S3
phi1 <- acos(runif(5000, -1,1))
phi2 <- runif(5000, 0,2*pi)
FigS3_dat = as.data.frame(do.call(rbind, lapply(1:5000, function(i){
  d <- cos(phi1[i])*t1+
    sin(phi1[i])*cos(phi2[i])*t2 +
    sin(phi1[i])*sin(phi2[i])*t3
  c(phi1[i],phi2[i],cor(d, t1), pcor(cbind(d, t2, t1))$estimate[1,2], pcor(cbind(d, t3, t1))$estimate[1,2])
})))
colnames(FigS3_dat) <- c("phi1","phi2", "corT1y", "corT2yT1","corT3yT1")
write.csv(FigS3_dat, file = "SourceData/FigureS3.csv", row.names = F)

jpeg("figures/FigureS3.jpeg", width = 21, height = 7, units = 'in', res = 600)
par(mfcol = c(1,3), oma = c(2,2,2,2), cex.axis = 0.8,cex = 2, cex.lab = 0.9, oma = c(0,0,0,0), cex.main = 1)
plot(FigS3_dat$phi1, FigS3_dat$corT1y,xlab = expression(paste(phi[1])), ylab = "Correlation coefficient", main = 
       expression('Correlation between t'[1]*' and y'))
mtext( "a", 2, adj=6, las=1, padj=-9, line =  -1, cex = 2.4, font = 2)
plot(FigS3_dat$phi2, FigS3_dat$corT2yT1,xlab = expression(paste(phi[2])), ylab = "Partial correlation coefficient", main = 
       expression('Partial correlation between t'[2]*' and y'))
mtext( "b", 2, adj=6, las=1, padj=-9, line =  -1, cex = 2.4, font = 2)
plot(FigS3_dat$phi2, FigS3_dat$corT3yT1,xlab = expression(paste(phi[2])), ylab = "Partial correlation coefficient", main = 
       expression('Partial correlation between t'[3]*' and y'))
  
mtext( "c", 2, adj=6, las=1, padj=-9, line =  -1, cex = 2.4, font = 2)
dev.off()

for (i in 1:3){
  plot(phi1, out[,1])
  mtext(letters[i], 2, adj=6, las=1, padj=-9, line =  -1, cex = 2.4, font = 2)
}
dev.off()


#### Figure S4
phi1 <- runif(5000, 0,pi)
phi2 <- runif(5000, 0,2*pi)
t1 <- rnorm(200)
t2 <- rnorm(200)
t3 <- rnorm(200)
out = do.call(rbind, lapply(1:2000, function(i){
  d <- cos(phi1[i])*t1+
    sin(phi1[i])*cos(phi2[i])*t2 +
    sin(phi1[i])*sin(phi2[i])*t3
  c(cor(d, t1), cor(d, t2), cor(d,t3))
}))

phi1 <- acos(runif(5000, -1,1))
phi2 <- runif(5000, 0,2*pi)
out2 = do.call(rbind, lapply(1:2000, function(i){
  d <- cos(phi1[i])*t1+
    sin(phi1[i])*cos(phi2[i])*t2 +
    sin(phi1[i])*sin(phi2[i])*t3
  c(cor(d, t1), cor(d, t2), cor(d,t3))
}))


jpeg("figures/FigureS4.jpeg", width = 21, height = 14, units = 'in', res = 600)
par(mfcol = c(2,3), oma = c(2,2,2,2), cex = 2, cex.lab = 1, oma = c(0,0,0,0), cex.main = 1.2)
for (i in 1:3){
  hist(out[,i], breaks = 50, freq = F,main = bquote("Correlation between t"[.(i)]*" and y"), xlab = "Correlation value")
  mtext(letters[i], 2, adj=6, las=1, padj=-9, line =  -1, cex = 2.4, font = 2)
  hist(out2[,i], breaks = 50, freq = F,main = bquote("Correlation between t"[.(i)]*" and y"), xlab = "Correlation value")
  mtext(letters[i+3], 2, adj=6, las=1, padj=-9, line =  -1, cex = 2.4, font = 2)
  
}
dev.off()
FigS4_dat <- data.frame(out, out2)
colnames(FigS4_dat) <- c("corT1_y_prior1", "corT2_y_prior1", "corT3_y_prior1", "corT1_y_prior2", "corT2_y_prior2", "corT3_y_prior2")
write.csv(FigS4_dat, file = "SourceData/FigureS4.csv", row.names = F)

### Figure S5
source("main/prior.R")
postDis <- read.csv(file = "SourceData/posterior.csv")
postDis <- read.csv(file = "SourceData/prior.csv")
priorDis <- sampler(nrow(postDis))
colnames(priorDis) <- list_params
write.csv(priorDis, file = "SourceData/prior.csv", row.names = F)

library(extraDistr)
titles <- c( expression(phi[paste(1,',', theta, 'min')]), expression(phi[paste(2,',', theta, 'min')]), expression(a[paste(theta,'min')]), expression(b[paste(theta,'min')]),
             expression(phi[paste(1,',l')]), expression(phi[paste(2,',l')]),expression(paste(a['l'])), expression(paste(b['l'])),
             expression(phi[paste(1,',c')]), expression(phi[paste(2,',c')]),expression(paste(a['c']))
)
names(titles) <- c("Tmin_phi1","Tmin_phi2", "Tmin_a","Tmin_b",
                   "l_phi1", "l_phi2", "l_a", "l_b",
                   "c_phi1","c_phi2", "c_a")

jpeg("figures/FigureS5.jpeg", width = 9, height = 12, units = "cm", res = 600)
par(mfcol= c(4,3), cex = 0.5, mar = c(4,2,1,1), oma = c(1,0.5,1,0)*2, cex.axis = 0.8, cex.lab = 1.2)
cex.side <- 0.65
cex.side2<- 0.5
for (i in c("Tmin_phi1","Tmin_phi2", "Tmin_a","Tmin_b",
            "l_phi1", "l_phi2", "l_a", "l_b",
            "c_phi1","c_phi2", "c_a")){
  dens.post <- stats::density(postDis[,i])
  dens.prior <- stats::density(priorDis[,i])
  plot(dens.post$x, dens.post$y, xlim = as.numeric(bounds[i,1:2]), main ="", xlab = titles[i], ylab ='', yaxt = "n", type = "l", col = "#CC6677", lwd = 2)
  points(dens.prior$x, dens.prior$y, type ="l", col = "#6699CC", lwd = 2)
  abline(v= bounds[i,], lty = 2)
  
  if (i == "Tmin_phi1"){
    mtext(side = 3, expression(paste("min. tol. temperature (", theta['min'], ')')), line = 1, cex = cex.side2)
    mtext(side = 2, expression(paste(phi[1])), line = 1, cex = cex.side)
  }
  if (i == "Tmin_phi2") mtext(side = 2, expression(paste(phi[2])), line = 1, cex = cex.side)
  if (i == "Tmin_a") mtext(side = 2, 'a', line = 1, cex = cex.side)
  if (i == "Tmin_b") mtext(side = 2, 'b', line = 1, cex = cex.side)
  
  if (i == "l_phi1") mtext(side = 3, "Sensitivity to biomass (l)", line = 1, cex = cex.side2)
  if (i == "c_phi1") mtext(side = 3, "Intra. comp. (c)", line = 1, cex = cex.side2)
}
plot(c(0,10), c(0,10),type  = "n",xaxt = "n", yaxt = "n", bty = "n", xlab ="", ylab = "n")
legend("topleft", legend = c("Prior", "Posterior"), cex =1.2, bty = "n", lwd = 2,lty = 1, col = c("#6699CC","#CC6677"), border = "white")
dev.off()
