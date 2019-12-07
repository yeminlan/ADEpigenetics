library(corrplot)
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)

mtl <- read.table('multiMTL.bed',sep="\t",header=F,stringsAsFactors = F)
mtl$ID <- paste0("MTL_",rownames(mtl))
colnames(mtl) <- c("chr","start","end","ID")

in.mtl <- read.table('InMTL.txt',sep="\t",header=T,stringsAsFactors = F)
mtl <- cbind(mtl,in.mtl)
rm(in.mtl)

auc <- read.table('AUC.txt',sep="\t",header=T,stringsAsFactors = F)
rownames(auc) <- mtl$ID

sample.list <- as.data.frame(colnames(auc))
colnames(sample.list) <- "sample"
sample.list$group <- as.factor(gsub("\\..*$","",sample.list$sample))
sample.list$marker <- as.factor(gsub("^.*\\.","",sample.list$sample))

mtl$CTCF.Y.auc <- log2( rowMeans(auc[,colnames(auc) %like% "Y.*.CTCF"]) + 1)
mtl$CTCF.O.auc <- log2( rowMeans(auc[,colnames(auc) %like% "O.*.CTCF"]) + 1)
mtl$CTCF.A.auc <- log2( rowMeans(auc[,colnames(auc) %like% "AD.*.CTCF"]) + 1)
mtl$Rad21.Y.auc <- log2( rowMeans(auc[,colnames(auc) %like% "Y.*.Rad21"]) + 1)
mtl$Rad21.O.auc <- log2( rowMeans(auc[,colnames(auc) %like% "O.*.Rad21"]) + 1)
mtl$Rad21.A.auc <- log2( rowMeans(auc[,colnames(auc) %like% "AD.*.Rad21"]) + 1)
mtl$H3K27ac.Y.auc <- log2( rowMeans(auc[,colnames(auc) %like% "Y.*.H3K27ac"]) + 1)
mtl$H3K27ac.O.auc <- log2( rowMeans(auc[,colnames(auc) %like% "O.*.H3K27ac"]) + 1)
mtl$H3K27ac.A.auc <- log2( rowMeans(auc[,colnames(auc) %like% "AD.*.H3K27ac"]) + 1)
mtl$H3K9ac.Y.auc <- log2( rowMeans(auc[,colnames(auc) %like% "Y.*.H3K9ac"]) + 1)
mtl$H3K9ac.O.auc <- log2( rowMeans(auc[,colnames(auc) %like% "O.*.H3K9ac"]) + 1)
mtl$H3K9ac.A.auc <- log2( rowMeans(auc[,colnames(auc) %like% "AD.*.H3K9ac"]) + 1)
mtl$H3K122ac.Y.auc <- log2( rowMeans(auc[,colnames(auc) %like% "Y.*.H3K122ac"]) + 1)
mtl$H3K122ac.O.auc <- log2( rowMeans(auc[,colnames(auc) %like% "O.*.H3K122ac"]) + 1)
mtl$H3K122ac.A.auc <- log2( rowMeans(auc[,colnames(auc) %like% "AD.*.H3K122ac"]) + 1)
mtl$H3K4me1.Y.auc <- log2( rowMeans(auc[,colnames(auc) %like% "Y.*.H3K4me1"]) + 1)
mtl$H3K4me1.O.auc <- log2( rowMeans(auc[,colnames(auc) %like% "O.*.H3K4me1"]) + 1)
mtl$H3K4me1.A.auc <- log2( rowMeans(auc[,colnames(auc) %like% "AD.*.H3K4me1"]) + 1)

write.csv(mtl,"MTL.result.csv",row.names = F)

