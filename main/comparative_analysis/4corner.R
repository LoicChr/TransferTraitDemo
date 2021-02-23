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
tab <- summary(Fourthcorner)[,c("Test", "Obs","Std.Obs","Alter","Pvalue","Pvalue.adj")]
tab[,1] <- gsub("V1 / ","", tab[,1])
tab[,2] <- gsub(" ","", tab[,2])
tab[,3] <- gsub(" ","", tab[,3])
write.table(tab, file = "SourceData/TabS6.csv", quote =F, row.names = F)
