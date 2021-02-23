##################################
#                                #
#     Fourth corner analysis     #
#                                #
##################################
setwd("./TransferTraitDemo/")

library(ade4)

### Dataset
load("data/data.Rdata")
env <- scale(temp)
env = as.data.frame(env)
# Main trait axes
spxt_log <- spxt
spxt_log[,1:6] <-apply(spxt_log[,1:6], 2, log)
spxt_log <- apply(spxt_log, 2, function(x){
  x[is.na(x)] <- mean(x, na.rm = T)
  x
})


Fourthcorner <- fourthcorner(as.data.frame(env), as.data.frame(ixp), as.data.frame(spxt_log))
tab <- summary(Fourthcorner)[,c("Test", "Obs","Pvalue","Pvalue.adj")]
tab[,"Test"] <- gsub("V1 / ","", tab[,"Test"])
tab[,"Obs"] <- gsub(" ","", tab[,"Obs"])
tab[,"Pvalue"] <- gsub(" ","", tab[,"Pvalue"])
tab[,1] <- c("Plant reproductive height (log)","Plant vegetative height (log)", "SLA  (log)", "LDMC (log)", "Leaf carbon content (log)", "Leaf nitrogen content (log)","Leaf d13C","Leaf d15N")

write.table(tab, file = "SourceData/TabS6.csv", quote =F, row.names = F)
