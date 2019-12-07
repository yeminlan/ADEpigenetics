d <- read.csv('MTL.result.log2.csv',check.names=F)


tmp <- read.table('EntorhinalCortex.sig.txt',header=T)
d$entor <- tmp[,1]
tmp <- read.table('PrefrontalCortex.sig.txt',header=F)
d$prefr <- tmp[,1]
tmp <- read.table('PrefrontalCortex_TauBurden.sig.txt',header=F)
d$prefr.tau <- tmp[,1]
rm(tmp)

tmp <- subset(d,select=c("H3K122ac.category","prefr.tau"))
tmp[,3] <- "others"
tmp[tmp[,1]=="DD_Gain" ,3] <- "DD_Gain"
tmp[tmp[,1]=="DD_Loss" ,3] <- "DD_Loss"
as.data.frame(table(paste0(tmp[,3],":",tmp[,2])))



write.table(subset(d,H3K27ac.category=="DD_Gain" & entor=="gain",select="Locus"),'~/tmp/DD_gain_EntorGain.bed',row.names=F,col.names=F,quote=F)

write.table(subset(d,H3K27ac.category=="DD_Loss" & entor=="loss",select="Locus"),'~/tmp/DD_loss_EntorLoss.bed',row.names=F,col.names=F,quote=F)

write.table(subset(d,H3K9ac.category=="DD_Gain" & prefr.tau=="gain",select="Locus"),'~/tmp/DD_gain_PrefrTauPos.bed',row.names=F,col.names=F,quote=F)

write.table(subset(d,H3K9ac.category=="DD_Loss" & prefr.tau=="loss",select="Locus"),'~/tmp/DD_loss_PrefrTauNeg.bed',row.names=F,col.names=F,quote=F)

