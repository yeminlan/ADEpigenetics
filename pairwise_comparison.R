library(corrplot)
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)

MTL <- read.csv('MTL.result.csv')
AUC <- read.table('AUC.txt',sep="\t",header=T,stringsAsFactors = F)

################ H3K27ac ################

auc <- AUC[,colnames(AUC) %like% "H3K27ac"]
mtl <- subset(MTL,select="ID")

sample.list <- as.data.frame(colnames(auc))
colnames(sample.list) <- "sample"
sample.list$group <- as.factor(gsub("\\..*$","",sample.list$sample))

## pairwise comparison of Y_vs_O
s1 <- sample.list[sample.list$group == "Y",]
s2 <- sample.list[sample.list$group == "O",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$Y.vs.O.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$Y.vs.O.qval <- p.adjust(mtl$Y.vs.O.pval,method="fdr")
mtl$Y.vs.O.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

## pairwise comparison of Y_vs_AD
s1 <- sample.list[sample.list$group == "Y",]
s2 <- sample.list[sample.list$group == "AD",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$Y.vs.AD.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$Y.vs.AD.qval <- p.adjust(mtl$Y.vs.AD.pval,method="fdr")
mtl$Y.vs.AD.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

## pairwise comparison of O_vs_AD
s1 <- sample.list[sample.list$group == "O",]
s2 <- sample.list[sample.list$group == "AD",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$O.vs.AD.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$O.vs.AD.qval <- p.adjust(mtl$O.vs.AD.pval,method="fdr")
mtl$O.vs.AD.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

mtl$ID <- NULL
colnames(mtl) <- paste0("H3K27ac.",colnames(mtl))
MTL <- cbind(MTL,mtl)

################ H3K9ac ################

auc <- AUC[,colnames(AUC) %like% "H3K9ac"]
mtl <- subset(MTL,select="ID")

sample.list <- as.data.frame(colnames(auc))
colnames(sample.list) <- "sample"
sample.list$group <- as.factor(gsub("\\..*$","",sample.list$sample))

## pairwise comparison of Y_vs_O
s1 <- sample.list[sample.list$group == "Y",]
s2 <- sample.list[sample.list$group == "O",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$Y.vs.O.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$Y.vs.O.qval <- p.adjust(mtl$Y.vs.O.pval,method="fdr")
mtl$Y.vs.O.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

## pairwise comparison of Y_vs_AD
s1 <- sample.list[sample.list$group == "Y",]
s2 <- sample.list[sample.list$group == "AD",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$Y.vs.AD.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$Y.vs.AD.qval <- p.adjust(mtl$Y.vs.AD.pval,method="fdr")
mtl$Y.vs.AD.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

## pairwise comparison of O_vs_AD
s1 <- sample.list[sample.list$group == "O",]
s2 <- sample.list[sample.list$group == "AD",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$O.vs.AD.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$O.vs.AD.qval <- p.adjust(mtl$O.vs.AD.pval,method="fdr")
mtl$O.vs.AD.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

mtl$ID <- NULL
colnames(mtl) <- paste0("H3K9ac.",colnames(mtl))
MTL <- cbind(MTL,mtl)

################ H3K122ac ################

auc <- AUC[,colnames(AUC) %like% "H3K122ac"]
mtl <- subset(MTL,select="ID")

sample.list <- as.data.frame(colnames(auc))
colnames(sample.list) <- "sample"
sample.list$group <- as.factor(gsub("\\..*$","",sample.list$sample))

## pairwise comparison of Y_vs_O
s1 <- sample.list[sample.list$group == "Y",]
s2 <- sample.list[sample.list$group == "O",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$Y.vs.O.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$Y.vs.O.qval <- p.adjust(mtl$Y.vs.O.pval,method="fdr")
mtl$Y.vs.O.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

## pairwise comparison of Y_vs_AD
s1 <- sample.list[sample.list$group == "Y",]
s2 <- sample.list[sample.list$group == "AD",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$Y.vs.AD.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$Y.vs.AD.qval <- p.adjust(mtl$Y.vs.AD.pval,method="fdr")
mtl$Y.vs.AD.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

## pairwise comparison of O_vs_AD
s1 <- sample.list[sample.list$group == "O",]
s2 <- sample.list[sample.list$group == "AD",]
x1 <- subset(auc,select=s1$sample)
x2 <- subset(auc,select=s2$sample)
## wilcox
for (i in 1:dim(x1)[1]) {
  mtl$O.vs.AD.pval[i] <- wilcox.test( as.numeric(x1[i,]), as.numeric(x2[i,]) )$p.value
}
mtl$O.vs.AD.qval <- p.adjust(mtl$O.vs.AD.pval,method="fdr")
mtl$O.vs.AD.diff <- rowMeans(x2) - rowMeans(x1)
rm(s1,s2,x1,x2,i)

mtl$ID <- NULL
colnames(mtl) <- paste0("H3K122ac.",colnames(mtl))
MTL <- cbind(MTL,mtl)

################################

write.csv(MTL,'MTL.result.csv',row.names = F)
