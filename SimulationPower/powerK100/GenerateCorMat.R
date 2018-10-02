## Generate K=20,40, 100 correlation matrices with specified eigenvalues 
rm(list=ls())
ls()
source("./GenCorr.R")
library(corrplot)

GennCorrLiu = function(K=20){
tmp = seq(200,1,length=K)
tmp[1] = 3000
eigenval1 = tmp/sum(tmp)*K
print(eigenval1)
print(sum(eigenval1))
mat = GenCorr(eigenval1,seed=123)
#fileName = paste0("CorMatK",as.character(K),".pdf")
#tl.cex = 1
#if(K==100){
#tl.cex = 0.5
#}
#pdf(fileName)

#corrplot(mat,tl.cex=tl.cex)
#dev.off()
print(eigen(mat)$values)
mat
}

###
sparse signals
###
K1=10
K0=90
K= K1+ K0
sigmaK0 = GennCorrLiu(K=K0)
sigmaK1 = GennCorrLiu(K=K1)
eigen(sigmaK1)$values
eigen(sigmaK0)$values


btw = matrix(runif(K1*K0,min=0,max = 0.001),K1,K0)
dim(btw)
sigma_part1 = cbind(sigmaK1,btw)
dim(sigma_part1)
dim(sigma_part1)
sigma_part2 = cbind(t(btw),sigmaK0)
dim(sigma_part2)
sigma = rbind(sigma_part1,sigma_part2)
dim(sigma)
eigen(sigma)$values
eigen(sigmaK1)$values
round(eigen(sigma)$vectors[,100],2)
eigen(sigmaK1)$vectors[,10]
sum(round(sigma,4)==t(round(sigma,4)))

write.csv(sigma, file=paste0("CorMatK",as.character(K),".csv"),row.names = F)
library(corrplot)
corrplot(round(sigma,4))
max(sigma)
min(sigma)
