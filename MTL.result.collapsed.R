## add CTCF/Rad21 within 5kb
d <- read.csv('MTL.result.log2.csv',header=T,check.names=F)
t1 <- read.table('MTL.hasCTCF_5kb.bed',header=F)
t2 <- read.table('MTL.hasRad21_5kb.bed',header=F)
d$hasCTCF.5kb <- t1$V1
d$hasRad21.5kb <- t2$V1
write.csv(d,'MTL.result.log2.csv',row.names=F)	

## collapse MTL table so that each transcript appear once (record only the nearest peak to its TSS)
d <- d[order(d$Distance,decreasing=F),]
d2 <- subset(d, !duplicated(Genes)) 
d2 <- d2[order(as.numeric(rownames(d2)),decreasing=F),]
write.csv(d2,'MTL.result.log2.collapsed.csv',row.names=F)

