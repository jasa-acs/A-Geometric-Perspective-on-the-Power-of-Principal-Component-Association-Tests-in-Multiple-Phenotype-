
###########################
## before LD pruning 
###########################

rm(list=ls())
hits = read.csv("./tophits/M8Union_newhits_PCO.csv",header = T)
sig = 5e-08
head(hits)
sig = 5e-08
hits$PC1 = as.numeric(hits$p.PC1 <sig )
hits$PC2 = as.numeric(hits$p.PC2 <sig )
hits$PC3 = as.numeric(hits$p.PC3 <sig )
hits$PC4 = as.numeric(hits$p.PC4 <sig )
hits$PC5 = as.numeric(hits$p.PC5 <sig )
hits$PC6 = as.numeric(hits$p.PC6 <sig )
hits$PC7 = as.numeric(hits$p.PC7 <sig )
hits$PC8 = as.numeric(hits$p.PC8 <sig )

#hits$PC10 = as.numeric(hits$p.PC10 <sig )

hits$PCMinP = as.numeric(hits$p.PCMinP <sig )
hits$PCFisher = as.numeric(hits$p.PCFisher <sig )
hits$PCLC = as.numeric(hits$p.PCLC <sig )
hits$WI = as.numeric(hits$p.WI <sig )
hits$Wald = as.numeric(hits$p.Wald <sig )
hits$VC = as.numeric(hits$p.VC <sig )
hits$PCAQ = as.numeric(hits$p.PCAQ <sig )
hits$PCO = as.numeric(hits$p.PCO <sig )
hits$tates = as.numeric(hits$p.tates <sig )
attach(hits)

#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)

membership1 = cbind(PC1,PC2,PC3)
membership2 = cbind(PC4,PC5,PC6)
membership3 = cbind(PC7,PC8,Wald)
membership4 = cbind(Wald,WI,VC)
membership5 = cbind(PC8,Wald,VC)
membership6 = cbind(PC1,Wald,VC)
membership7 = cbind(PCLC,Wald,PCFisher)
membership8 = cbind(Wald,VC,PCAQ)
membership9 = cbind(Wald,PCLC,PCO)
membership10 = cbind(PCAQ,PCMinP,PCO)

a1 = vennCounts(membership1)
a2 = vennCounts(membership2)
a3 = vennCounts(membership3)
a4 = vennCounts(membership4)
a5 = vennCounts(membership5)
a6 = vennCounts(membership6)
a7 = vennCounts(membership7)
a8 = vennCounts(membership8)
a9 = vennCounts(membership9)
a10 = vennCounts(membership10)

pdf("Venn_M8_Not_pruned.pdf",height = 8,width = 18)
par(mfrow=c(2,5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
par(mar = c(0, 0, 0, 0), oma = c(0.5, 0.5, 0.5, 0.5))
vennDiagram(a1,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a2,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a3,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a4,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a5,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a6,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a7,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a8,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a9,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
vennDiagram(a10,cex=1.2,mar = c(0, 0, 0, 0),circle.col=1,counts.col=1)
dev.off()

#################
##### counting 
################
dat = hits
ind.PC1 = which(dat$p.PC1 < sig)
ind.PC2 = which(dat$p.PC2 < sig)
ind.PC3 = which(dat$p.PC3 < sig)
ind.PC4 = which(dat$p.PC4 < sig)

ind.PC5 = which(dat$p.PC5 < sig)
ind.PC6 = which(dat$p.PC6 < sig)
ind.PC7 = which(dat$p.PC7 < sig)
ind.PC8 = which(dat$p.PC8 < sig)



ind.PCMinP = which(dat$p.PCMinP < sig)
ind.PCFisher = which(dat$p.PCFisher < sig)
ind.PCLC = which(dat$p.PCLC < sig)
ind.WI = which(dat$p.WI < sig)
ind.Wald = which(dat$p.Wald < sig)
ind.VC = which(dat$p.VC < sig)

ind.PCAQ = which(dat$p.PCAQ < sig)
ind.PCO = which(dat$p.PCO < sig)
ind.tates = which(dat$p.tates < sig)

ind.list = list(ind.PC1,ind.PC2,ind.PC3,ind.PC4,ind.PC5,ind.PC6,ind.PC7,ind.PC8,ind.PCMinP,ind.PCFisher,ind.PCLC,ind.WI,ind.Wald,ind.VC,ind.PCAQ,ind.PCO,ind.tates)
num = unlist(lapply(ind.list,length))
names(num)=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PCMinP","PCFisher","PCLC","WI","Wald","VC","PCAQ","PCO","tates")
write.csv(t(num),file="M8_Not_pruned_count.csv",row.names=F)