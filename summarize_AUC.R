library(plyr)
library(data.table)

## read AUC.txt from single mark MTL
d <- read.table('t3',header=F,sep="\t")
# multiMTL
dA <- d[,1:3]
colnames(dA) <- c("chr","start","end")
# MTL
dB <- d[,4:6]
colnames(dB) <- c("chr","start","end")
# AUC
dC <- d[,7:dim(d)[2]]
h <- read.table('t2',header=T,sep="\t",comment.char="",check.names=F)
colnames(dC) <- colnames(h)[4:dim(h)[2]]
rm(h)

# use AUC*width to aggregate
dA$coordinate <- paste0(dA$chr,":",dA$start,"-",dA$end)
dB$width <- dB$end-dB$start
dC <- dC * dB$width/1000
d <- cbind( subset(dA,select="coordinate") , subset(dB,select="width") , dC )
D <- aggregate(d[,-1] ,by=list(d$coordinate),FUN="sum")
colnames(D)[1] <- "coordinate"

## divide aggregated sum by aggregated width
DA <- D[,1:2]
DA$chr <- gsub(":.*$","",DA$coordinate)
DA$start <- gsub("-.*$","",gsub(".*:","",DA$coordinate))
DA$end <- gsub(".*-","",DA$coordinate)
DB <- D[,3:dim(D)[2]]
DB <- DB*1000/DA$width
D <- cbind( subset(DA,select=c("coordinate")) ,DB)
# 

## join to multiMTL table and replace zeros
mtl <- read.table('multiMTL.bed',header=F,sep="\t")
colnames(mtl) <- c("chr","start","end")
mtl$coordinate <- paste0(mtl$chr,":",mtl$start,"-",mtl$end)
mtl <- join(mtl,D)
D <- mtl[,c(-1:-4)]
D[is.na(D)] <- 0

write.table(D,'AUC.txt',sep="\t",row.names=F,quote=F)


