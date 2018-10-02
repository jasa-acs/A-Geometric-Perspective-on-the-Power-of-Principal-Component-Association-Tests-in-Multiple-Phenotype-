############################################################
## PCA demo for significant and not-significant points    ##
############################################################
rm(list=ls())
library(mvtnorm)
#library(MPAT)
source("./tatespvalZ.R")
tatespvalZ()

genSample = function(n,mu,Sigma){
	Z = rmvnorm(n,mean=mu,sigma=Sigma)
	return(Z)	
}

addPoint = function(){
points(0,0,pch="+",col="red",cex=2)
}

addLines = function(){
abline(h=1.96,col="blue",lty=3)
abline(v=1.96,col="blue",lty=3)
abline(h=-1.96,col="blue",lty=3)
abline(v=-1.96,col="blue",lty=3)
}
addArrows = function(){

if(rho > 0) {
arrows(x0=-3*sqrt(lambda1),y=-3*sqrt(lambda1),x1=3*sqrt(lambda1),y1=3*sqrt(lambda1),code=3,col="blue",length=0.15)
arrows(x0=3*sqrt(lambda2),y=-3*sqrt(lambda2),x1=-3*sqrt(lambda2),y1=3*sqrt(lambda2),code=3,col="green",length=0.15)
}
if(rho < 0){
arrows(x0=-3*sqrt(lambda2),y=-3*sqrt(lambda2),x1=3*sqrt(lambda2),y1=3*sqrt(lambda2),code=3,col="blue",length=0.15)
arrows(x0=3*sqrt(lambda1),y=-3*sqrt(lambda1),x1=-3*sqrt(lambda1),y1=3*sqrt(lambda1),code=3,col="green",length=0.15)

}

}




############################
############################
rho = 0.6
Sigma=cbind(c(1,rho),c(rho,1))
Sigma
u=eigen(Sigma)$vectors
lambda = eigen(Sigma)$values
lambda1 = lambda[1]
lambda2 = lambda[2]
J = as.matrix(c(1,1))

outfilename = paste0("TATES_RejectionBoundary.png")
png(outfilename)
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = tatespvalZ(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="red",cex.lab=1,lwd=1,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "TATES",cex.main=2,font.main=1)
addPoint()
addLines()
addArrows()
dev.off()
