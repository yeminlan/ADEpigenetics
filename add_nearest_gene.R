library(plyr)
library(data.table)

mtl <- read.csv('MTL.result.csv',stringsAsFactors = F)

t <- read.table('multiMTL.anno.txt',header=T,sep="\t",comment.char = "",quote = "",stringsAsFactors = F)
t <- t[,c(2,3,4,10,16)]
colnames(t) <- c("chr","start","end","dist.to.nearest.gene","nearest.gene")
t$start <- t$start-1
mtl <- join(mtl,t)
mtl <- mtl[,c( 1:3,dim(mtl)[2]-1,dim(mtl)[2],4:(dim(mtl)[2]-2) )]

## merge coordinate 3-columns to 1 column
mtl$coordinate <- paste0(mtl$chr,":",mtl$start,"-",mtl$end)
mtl <- mtl[, c( dim(mtl)[2],4:(dim(mtl)[2]-1) ) ]

##
write.csv(mtl,'MTL.result.csv',row.names = F)
