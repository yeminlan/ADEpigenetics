mtl <- read.table('t1',header=F)
colnames(mtl) <- 'Distance'
mtl$Distance.group <- cut(mtl$Distance,breaks=c(0,1000,50000,100000,Inf),include.lowest=T)
levels(mtl$Distance.group) <- c("within.1kb","1-50kb","50-100kb","beyond.100kb")
write.csv(mtl,'distance_to_ensembleTSS.csv',row.names=F,quote=F)


