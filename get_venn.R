library(plyr)
library(data.table)
library(venneuler)

mtl <- read.csv('MTL.result.csv',stringsAsFactors = F)

###########################
## replace with ensembleTSS
#d <- read.csv('distance_to_ensembleTSS.csv',header=T)
#mtl$dist.to.nearest.gene <- d$Distance
#rm(d)
###########################

mtl$group1 <- "within.1kb"
mtl$group1[abs(mtl$dist.to.nearest.gene)>1000] <- ">1kb"
mtl$group2 <- "no.H3K4me1"
mtl$group2[(mtl$In.H3K4me1.Y+mtl$In.H3K4me1.O+mtl$In.H3K4me1.A)>0] <- "has.H3K4me1"
table(subset(mtl,select=c("group1","group2")))

#################

mtl.part <- subset(mtl,group1=="within.1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t1 <- as.data.frame(table(mtl.part$group))
colnames(t1) <- c("group","Y")
mtl.part <- subset(mtl,group1=="within.1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t2 <- as.data.frame(table(mtl.part$group))
colnames(t2) <- c("group","O")
mtl.part <- subset(mtl,group1=="within.1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t3 <- as.data.frame(table(mtl.part$group))
colnames(t3) <- c("group","A")
mtl.part <- subset(mtl,group1=="within.1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y","In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O","In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$In.H3K27ac <- mtl.part$In.H3K27ac.Y + mtl.part$In.H3K27ac.O + mtl.part$In.H3K27ac.A
mtl.part$In.H3K9ac <- mtl.part$In.H3K9ac.Y + mtl.part$In.H3K9ac.O + mtl.part$In.H3K9ac.A
mtl.part$In.H3K122ac <- mtl.part$In.H3K122ac.Y + mtl.part$In.H3K122ac.O + mtl.part$In.H3K122ac.A
mtl.part <- subset(mtl.part,select=c("In.H3K27ac","In.H3K9ac","In.H3K122ac"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t4 <- as.data.frame(table(mtl.part$group))
colnames(t4) <- c("group","any")
t <- join(t1,t2)
t <- join(t,t3)
t <- join(t,t4)

#################

mtl.part <- subset(mtl,group1=="within.1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t1 <- as.data.frame(table(mtl.part$group))
colnames(t1) <- c("group","Y")
mtl.part <- subset(mtl,group1=="within.1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t2 <- as.data.frame(table(mtl.part$group))
colnames(t2) <- c("group","O")
mtl.part <- subset(mtl,group1=="within.1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t3 <- as.data.frame(table(mtl.part$group))
colnames(t3) <- c("group","A")
mtl.part <- subset(mtl,group1=="within.1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y","In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O","In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$In.H3K27ac <- mtl.part$In.H3K27ac.Y + mtl.part$In.H3K27ac.O + mtl.part$In.H3K27ac.A
mtl.part$In.H3K9ac <- mtl.part$In.H3K9ac.Y + mtl.part$In.H3K9ac.O + mtl.part$In.H3K9ac.A
mtl.part$In.H3K122ac <- mtl.part$In.H3K122ac.Y + mtl.part$In.H3K122ac.O + mtl.part$In.H3K122ac.A
mtl.part <- subset(mtl.part,select=c("In.H3K27ac","In.H3K9ac","In.H3K122ac"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t4 <- as.data.frame(table(mtl.part$group))
colnames(t4) <- c("group","any")
t <- join(t1,t2)
t <- join(t,t3)
t <- join(t,t4)

#################

mtl.part <- subset(mtl,group1==">1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t1 <- as.data.frame(table(mtl.part$group))
colnames(t1) <- c("group","Y")
mtl.part <- subset(mtl,group1==">1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t2 <- as.data.frame(table(mtl.part$group))
colnames(t2) <- c("group","O")
mtl.part <- subset(mtl,group1==">1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t3 <- as.data.frame(table(mtl.part$group))
colnames(t3) <- c("group","A")
mtl.part <- subset(mtl,group1==">1kb" & group2=="has.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y","In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O","In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$In.H3K27ac <- mtl.part$In.H3K27ac.Y + mtl.part$In.H3K27ac.O + mtl.part$In.H3K27ac.A
mtl.part$In.H3K9ac <- mtl.part$In.H3K9ac.Y + mtl.part$In.H3K9ac.O + mtl.part$In.H3K9ac.A
mtl.part$In.H3K122ac <- mtl.part$In.H3K122ac.Y + mtl.part$In.H3K122ac.O + mtl.part$In.H3K122ac.A
mtl.part <- subset(mtl.part,select=c("In.H3K27ac","In.H3K9ac","In.H3K122ac"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t4 <- as.data.frame(table(mtl.part$group))
colnames(t4) <- c("group","any")
t <- join(t1,t2)
t <- join(t,t3)
t <- join(t,t4)

#################

mtl.part <- subset(mtl,group1==">1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t1 <- as.data.frame(table(mtl.part$group))
colnames(t1) <- c("group","Y")
mtl.part <- subset(mtl,group1==">1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t2 <- as.data.frame(table(mtl.part$group))
colnames(t2) <- c("group","O")
mtl.part <- subset(mtl,group1==">1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t3 <- as.data.frame(table(mtl.part$group))
colnames(t3) <- c("group","A")
mtl.part <- subset(mtl,group1==">1kb" & group2=="no.H3K4me1",select=c("In.H3K27ac.Y","In.H3K9ac.Y","In.H3K122ac.Y","In.H3K27ac.O","In.H3K9ac.O","In.H3K122ac.O","In.H3K27ac.A","In.H3K9ac.A","In.H3K122ac.A"))
mtl.part$In.H3K27ac <- mtl.part$In.H3K27ac.Y + mtl.part$In.H3K27ac.O + mtl.part$In.H3K27ac.A
mtl.part$In.H3K9ac <- mtl.part$In.H3K9ac.Y + mtl.part$In.H3K9ac.O + mtl.part$In.H3K9ac.A
mtl.part$In.H3K122ac <- mtl.part$In.H3K122ac.Y + mtl.part$In.H3K122ac.O + mtl.part$In.H3K122ac.A
mtl.part <- subset(mtl.part,select=c("In.H3K27ac","In.H3K9ac","In.H3K122ac"))
mtl.part$group <- paste0(as.numeric(mtl.part[,1]>0),"-",as.numeric(mtl.part[,2]>0),"-",as.numeric(mtl.part[,3]>0))
t4 <- as.data.frame(table(mtl.part$group))
colnames(t4) <- c("group","any")
t <- join(t1,t2)
t <- join(t,t3)
t <- join(t,t4)

