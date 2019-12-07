library(plyr)
library(data.table)

mtl <- read.csv('MTL.result.csv', check.names=F, stringsAsFactors = F)
colnames(mtl) <- c("Locus","Distance","Genes","ID",
                   "CTCF.InA","CTCF.InO","CTCF.InY","Rad21.InA","Rad21.InO","Rad21.InY",
                   "H3K27ac.InA","H3K27ac.InO","H3K27ac.InY","H3K9ac.InA","H3K9ac.InO","H3K9ac.InY",
                   "H3K122ac.InA","H3K122ac.InO","H3K122ac.InY","H3K4me1.InA","H3K4me1.InO","H3K4me1.InY",
                   "CTCF.HY","CTCF.HO","CTCF.HA","Rad21.HY","Rad21.HO","Rad21.HA",
                   "H3K27ac.HY","H3K27ac.HO","H3K27ac.HA","H3K9ac.HY","H3K9ac.HO","H3K9ac.HA",
                   "H3K122ac.HY","H3K122ac.HO","H3K122ac.HA","H3K4me1.HY","H3K4me1.HO","H3K4me1.HA",
                   "H3K27ac.PHO-HY","H3K27ac.QHO-HY","H3K27ac.HO-HY","H3K27ac.PHA-HY","H3K27ac.QHA-HY","H3K27ac.HA-HY",
                   "H3K27ac.PHA-HO","H3K27ac.QHA-HO","H3K27ac.HA-HO","H3K9ac.PHO-HY","H3K9ac.QHO-HY","H3K9ac.HO-HY",
                   "H3K9ac.PHA-HY","H3K9ac.QHA-HY","H3K9ac.HA-HY","H3K9ac.PHA-HO","H3K9ac.QHA-HO","H3K9ac.HA-HO",
                   "H3K122ac.PHO-HY","H3K122ac.QHO-HY","H3K122ac.HO-HY","H3K122ac.PHA-HY","H3K122ac.QHA-HY","H3K122ac.HA-HY",
                   "H3K122ac.PHA-HO","H3K122ac.QHA-HO","H3K122ac.HA-HO")

mtl$Distance <- abs(mtl$Distance)
mtl$Distance.group <- cut(mtl$Distance,breaks=c(0,1000,50000,100000,Inf),include.lowest=T)
levels(mtl$Distance.group) <- c("within.1kb","1-50kb","50-100kb","beyond.100kb")

mtl$`H3K27ac.HO-HY` <- mtl$H3K27ac.HO - mtl$H3K27ac.HY
mtl$`H3K27ac.HA-HY` <- mtl$H3K27ac.HA - mtl$H3K27ac.HY
mtl$`H3K27ac.HA-HO` <- mtl$H3K27ac.HA - mtl$H3K27ac.HO
mtl$`H3K9ac.HO-HY` <- mtl$H3K9ac.HO - mtl$H3K9ac.HY
mtl$`H3K9ac.HA-HY` <- mtl$H3K9ac.HA - mtl$H3K9ac.HY
mtl$`H3K9ac.HA-HO` <- mtl$H3K9ac.HA - mtl$H3K9ac.HO
mtl$`H3K122ac.HO-HY` <- mtl$H3K122ac.HO - mtl$H3K122ac.HY
mtl$`H3K122ac.HA-HY` <- mtl$H3K122ac.HA - mtl$H3K122ac.HY
mtl$`H3K122ac.HA-HO` <- mtl$H3K122ac.HA - mtl$H3K122ac.HO

mtl <- subset(mtl,select=c("Locus","Genes","Distance","Distance.group","ID", colnames(mtl)[!(colnames(mtl) %in% c("Locus","Genes","Distance","Distance.group","ID"))] ))

rna <- read.csv('brain.rna.csv', check.names=F, stringsAsFactors = F)
colnames(rna)[1] <- "Genes"

mtl <- join(mtl,rna)

write.csv(mtl,'MTL.result.log2.csv',row.names = F)

