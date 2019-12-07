library(corrplot)
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)

auc <- read.table('AUC.txt',sep="\t",header=T,stringsAsFactors = F)
mtl <- read.table('MTL.txt',sep="\t",header=F,stringsAsFactors = F)
mtl$ID <- paste0("MTL_",rownames(mtl))
rownames(auc) <- mtl$ID
colnames(mtl)[1:6] <- c("chr","start","end","peak.in.AD","peak.in.O","peak.in.Y")

sample.list <- as.data.frame(colnames(auc))
colnames(sample.list) <- "sample"
sample.list$group <- as.factor(gsub("\\..*$","",sample.list$sample))

NeuFrac <- read.table('NeuFrac.txt',sep="\t",header=T,stringsAsFactors = F)
NeuFrac$sampleID <- gsub("-",".",NeuFrac$sampleID)

###############################

## spearman's correlation
pdf("correlation.pdf",height=8,width=7)
c <- cor(auc, method = "spearman")
col <- colorRampPalette(c("blue","grey90","red", "grey90", "blue"))
corrplot(c, type = "lower", method = "number", order="hclust", hclust.method = "ward.D", mar=c(0,3,0,0), 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.4, addshade = "all", 
         bg="white", cl.lim=c(0,1), col=col(100))
dev.off()
rm(c,col)

###############################

## mask peaks with 10% highest pearson corr with NeuFrac
c <- apply(auc, 1, function(x) cor(NeuFrac$percentage, x, method = "pearson") )
c <- abs(c)
mtl$mask <- as.numeric(c>quantile(c, probs = 0.9))

auc.filter <- auc[mtl$mask==0,]

## PCA plot of all/top10000 MTLs
pdf("PCA.mask.pdf",height=8,width=8)
#all MTLs
pc <- prcomp(t(auc.filter), scale=TRUE)
scores <- as.data.frame(pc$x)
v <- as.integer(100*(pc$sdev)^2/sum(pc$sdev^2))
for (i in 1:5){
  print( paste0("PC",i," (",v[i],"%): ",cor(NeuFrac$percentage, scores[,i], method = "pearson")) )
}
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), col=sample.list$group)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(alpha = 0.8, size = 2) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position = "none") +
  ggtitle("PCA plot of all MTLs") +
  xlab(paste0("PC1 (",v[1],"%)")) + 
  ylab(paste0("PC2 (",v[2],"%)"))
#top10000 MTLs
x <- auc.filter[order(rowSums(auc.filter),decreasing=T),]
x <- x[1:10000,]
pc <- prcomp(t(x), scale=TRUE)
scores <- as.data.frame(pc$x)
v <- as.integer(100*(pc$sdev)^2/sum(pc$sdev^2))
for (i in 1:5){
  print( paste0("PC",i," (",v[i],"%): ",cor(NeuFrac$percentage, scores[,i], method = "pearson")) )
}
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), col=sample.list$group)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(alpha = 0.8, size = 2) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position = "none") +
  ggtitle("PCA plot of top 10000 MTLs") +
  xlab(paste0("PC1 (",v[1],"%)")) + 
  ylab(paste0("PC2 (",v[2],"%)"))
dev.off()
rm(pc,scores,v,x)

###############################

## PCA plot of all/top10000 MTLs
pdf("PCA.pdf",height=8,width=8)
#all MTLs
pc <- prcomp(t(auc), scale=TRUE)
scores <- as.data.frame(pc$x)
v <- as.integer(100*(pc$sdev)^2/sum(pc$sdev^2))
for (i in 1:5){
  print( paste0("PC",i," (",v[i],"%): ",cor(NeuFrac$percentage, scores[,i], method = "spearman")) )
}
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), col=sample.list$group)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(alpha = 0.8, size = 2) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position = "none") +
  ggtitle("PCA plot of all MTLs") +
  xlab(paste0("PC1 (",v[1],"%)")) + 
  ylab(paste0("PC2 (",v[2],"%)"))
#top10000 MTLs
x <- auc[order(rowSums(auc),decreasing=T),]
x <- x[1:10000,]
pc <- prcomp(t(x), scale=TRUE)
scores <- as.data.frame(pc$x)
v <- as.integer(100*(pc$sdev)^2/sum(pc$sdev^2))
for (i in 1:5){
  print( paste0("PC",i," (",v[i],"%): ",cor(NeuFrac$percentage, scores[,i], method = "spearman")) )
}
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), col=sample.list$group)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(alpha = 0.8, size = 2) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position = "none") +
  ggtitle("PCA plot of top 10000 MTLs") +
  xlab(paste0("PC1 (",v[1],"%)")) + 
  ylab(paste0("PC2 (",v[2],"%)"))
dev.off()
rm(pc,scores,v,x)

###############################

## anova
for (i in 1:dim(auc)[1]) {
  d <- melt(auc[i,])
  colnames(d) <- c("sample","auc")
  d <- join(d,sample.list,by="sample")
  fit <- lm(auc ~ group, data=d)
  fit2 <- anova(fit)
  mtl$anova[i] <- fit2$`Pr(>F)`[1]
}
rm(d,i,fit,fit2)

## PCA plot of significant MTLs
auc.sig <- auc[mtl$anova<0.05,]
pdf("PCA.sig.pdf",height=8,width=8)
pc <- prcomp(t(auc.sig), scale=TRUE)
scores <- as.data.frame(pc$x)
v <- as.integer(100*(pc$sdev)^2/sum(pc$sdev^2))
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), col=sample.list$group)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(alpha = 0.8, size = 2) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position = "none") +
  ggtitle("PCA plot of significant MTLs") +
  xlab(paste0("PC1 (",v[1],"%)")) + 
  ylab(paste0("PC2 (",v[2],"%)"))
dev.off()
rm(pc,scores,v)

## k-means for anova.sig.MTLs
pdf("k-means_decide.pdf",height=4,width=6)
wss <- (nrow(auc.sig)-1)*sum(apply(auc.sig,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(auc.sig,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()

bestK <- 10

set.seed(1) 
wss <- kmeans(auc.sig,centers=bestK) 
d <- as.data.frame(wss$cluster)
d$ID <- row.names(d)
colnames(d) <- c("k.means","ID")
mtl <- join(mtl,d)
mtl$k.means[is.na(mtl$k.means)] <- 0 #non-sig MTLs were assigned to cluster0
rm(d,wss)

for (i in 1:bestK){
  t <- mtl[mtl$k.means==i,]
  print(c(wilcox.test(t$Y.auc,t$O.auc)$p.value, wilcox.test(t$Y.auc,t$AD.auc)$p.value, wilcox.test(t$O.auc,t$AD.auc)$p.value))
}
rm(i,t)

###############################

mtl$AD.auc <- rowMeans(auc[,colnames(auc) %like% "AD."])
mtl$O.auc <- rowMeans(auc[,colnames(auc) %like% "O."])
mtl$Y.auc <- rowMeans(auc[,colnames(auc) %like% "Y."])

###############################

d <- melt(mtl,id.vars=c("chr","start","end","peak.in.AD","peak.in.O","peak.in.Y","ID","anova","k.means"),
          variable.name = "group", value.name = "ave.AUC")
d <- d[d$k.means!=0,]
d$group <- gsub(".auc","",d$group)
d$group <- factor(d$group,c("Y","O","AD"))
pdf("k-means_boxplot.pdf",height=4,width=8)
ggplot(d,aes(group,ave.AUC,fill=group)) + geom_boxplot() + 
  facet_wrap(~k.means, nrow = 2, scales = "free") +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position = "none") + xlab("")
dev.off()
table(mtl$k.means)
rm(d)

###############################

write.csv(mtl,"MTL.result.csv",row.names = F)
