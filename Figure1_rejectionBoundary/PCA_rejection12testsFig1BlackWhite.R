############################################################
## PCA demo for significant and not-significant points    ##
############################################################
library(mvtnorm)
library(MPAT)
genSample = function(n,mu,Sigma){
	Z = rmvnorm(n,mean=mu,sigma=Sigma)
	return(Z)	
}

addPoint = function(){
points(0,0,pch="+",col="black",cex=2)
}

addLines = function(){
#abline(h=1.96,col="black",lty=3)
#abline(v=1.96,col="black",lty=3)
#abline(h=-1.96,col="black",lty=3)
#abline(v=-1.96,col="black",lty=3)
}
addArrows = function(){

if(rho > 0) {
arrows(x0=-3*sqrt(lambda1),y=-3*sqrt(lambda1),x1=3*sqrt(lambda1),y1=3*sqrt(lambda1),code=3,col="black",lty=1,lwd=2,length=0.15)
arrows(x0=3*sqrt(lambda2),y=-3*sqrt(lambda2),x1=-3*sqrt(lambda2),y1=3*sqrt(lambda2),code=3,col="black",lty=1,lwd=2,length=0.15)
}
if(rho < 0){
arrows(x0=-3*sqrt(lambda2),y=-3*sqrt(lambda2),x1=3*sqrt(lambda2),y1=3*sqrt(lambda2),code=3,col="black",lty=2,lwd=2,length=0.15)
arrows(x0=3*sqrt(lambda1),y=-3*sqrt(lambda1),x1=-3*sqrt(lambda1),y1=3*sqrt(lambda1),code=3,col="black",lty=3,lwd=2,length=0.15)

}

}


PCAQEmpDist = function(simNum=1e5,Sigma=Sigma){
# Step 1: simulation Z-vectors under the null
# equivalent to parametric Bootstrap
K = dim(Sigma)[1]
## this multivariate normal must be under the Null!!
## This is the basis for generating a null reference distribution
Z.mat = rmvnorm(simNum,rep(0,K),sigma = Sigma)
PCQMinP = rep(NA,simNum)
for(b in 1:simNum){
p.PCBall = WI(Z.mat[b,],Sigma)
p.PCSS = Wald(Z.mat[b,],Sigma)
p.PCWSS = VC(Z.mat[b,],Sigma)
## compute PCBall, PCSS, PCWSS, take minimum
PCQMinP[b] = min(c(p.PCBall,p.PCSS,p.PCWSS))
}
## return the empirical null distribution
return(PCQMinP)
}


PCMetaEmpDist = function(simNum=1e5,Sigma=Sigma){
# Step 1: simulation Z-vectors under the null
# equivalent to parametric Bootstrap
K = dim(Sigma)[1]
## this multivariate normal must be under the Null!!
## This is the basis for generating a null reference distribution
Z.mat = rmvnorm(simNum,rep(0,K),sigma = Sigma)
PCMetaMinP = rep(NA,simNum)
for(b in 1:simNum){
p.PCMinP = PCMinP(Z.mat[b,],Sigma)
p.PCFisher = PCFisher(Z.mat[b,],Sigma)
p.PCLC = PCLC(Z.mat[b,],Sigma)
p.WI = WI(Z.mat[b,],Sigma)
p.Wald = Wald(Z.mat[b,],Sigma)
p.VC = VC(Z.mat[b,],Sigma)
## compute PCBall, PCSS, PCWSS, take minimum
PCMetaMinP[b] = min(c(p.PCMinP,p.PCFisher,p.PCLC,p.WI,p.Wald,p.VC))
}
## return the empirical null distribution
return(PCMetaMinP)
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

outfilename = paste0("PCA_RejectionBoundary12testsBlackWhite.pdf")
pdf(outfilename,height=12,width=16.5)
#mat <- matrix(seq(1,12), nrow=3,byrow = TRUE)
#layout(mat)
par(mfrow=c(3,4))
par(mar=c(5.1,5.1,2.1,2.1))

#################
### the Oracle Test depends on alternative
#################
## suppose beta.alternative = c(2,0), a bit confusing 
###########################################
Oracle = function(Z.vec,Sigma,beta.alternative){
T = (t(beta.alternative)%*%solve(Sigma)%*%Z.vec)^2 /(t(beta.alternative)%*%solve(Sigma)%*%beta.alternative)
p = pchisq(T,df=1,lower.tail=F)
return(p)
}

PCMix = function(Z.vec,Sigma,method="davies"){

   U = eigen(Sigma)$vectors
   lambdas = eigen(Sigma)$values
   Z_star = t(U)%*%Z.vec
   Sigma_star = diag(lambdas)
   p = mixSD(Z_star,Sigma_star)
   return(p)
}

#####################
#### figure setup parameters
lwd=2
cex.lab=2
lty=2
cex.axis=2
font.lab=2
font.main=1
cex.main=2.5
#######################
# Oracle 1: beta = UJ, identical to PCLC 
##################
C = 0.05
n=100
x = seq(-10,10,length=n)
y = seq(-10,10,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = Oracle(c(x[i],y[j]),Sigma,beta.alternative=c(0,1))
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title(expression(paste("Oracle for ",beta," = (0,1)")),cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

########################
## Oracle 2: beta = (2,2), identical to PC1 
##################
C = 0.05
n=100
x = seq(-10,10,length=n)
y = seq(-10,10,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = Oracle(c(x[i],y[j]),Sigma,beta.alternative=c(1,1))
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title(expression(paste("Oracle for ",beta," = (1,1)")),cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

#################
### PCMix
#################
# C = 0.05
# n=100
# x = seq(-5,5,length=n)
# y = seq(-5,5,length=n)
# J = as.matrix(c(1,1))
# a = t(J)%*%solve(Sigma)

# z = matrix(NA,ncol=n,nrow=n)
# for(i in 1:n){
#    for(j in 1:n){
#     z[i,j] = PCMix(c(x[i],y[j]),Sigma)
#    }
# }
# contour(x,y,z,levels=C,drawlabels=FALSE,col="black",cex.lab=0.85,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
# title( "PCMix",cex.main=2,font.main=1)
# addPoint()
# addLines()
# addArrows()

#################
### PC1
#################
C = 0.05
n=100
x = seq(-10,10,length=n)
y = seq(-10,10,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = PC(c(x[i],y[j]),Sigma,1)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PC1",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

#################
### PC2
#################
C = 0.05
n=100
x = seq(-10,10,length=n)
y = seq(-10,10,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = PC(c(x[i],y[j]),Sigma,2)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PC2",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

#################
### PCMinP
#################
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = PCMinP(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PCMinP",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

#################
### PCFisher
#################
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = PCFisher(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PCFisher",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

#################
### PCLC
#################
C = 0.05
n=100
x = seq(-10,10,length=n)
y = seq(-10,10,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = PCLC(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PCLC",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()
w1 = 1/lambda1*1/sqrt(2) + 1/lambda2*(-1/sqrt(2))
w2 = 1/lambda1*1/sqrt(2) + 1/lambda2*(1/sqrt(2))
arrows(x0=2,y0=2*w2/w1,x1=-2,y1=-2*w2/w1,code=3,col="black",length=0.15,angle=90)
# arrows(x0=-2*w1/w2,y0=2,x1=-2*w1/w2,y1=-2,code=3,col="black",length=0.15)
#abline(0,-1/0.6,col="black")
#arrows(x0=2,y0=2*w2/w1,x1=-2,y1=-2*w2/w1,code=3,col="black",length=0.15)
#################
### WI
#################
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = WI(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "WI",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

#################
### Wald 
#################
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = Wald(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "Wald",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()


#################
### VC
#################
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
   for(j in 1:n){
    z[i,j] = VC(c(x[i],y[j]),Sigma)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "VC",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

##########
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)

z = matrix(NA,ncol=n,nrow=n)
PCAQ_ref = PCAQEmpDist(simNum=1e5, Sigma=Sigma)
for(i in 1:n){
   for(j in 1:n){
    z.WI = WI(c(x[i],y[j]),Sigma)
	z.Wald = Wald(c(x[i],y[j]),Sigma)
	z.VC= VC(c(x[i],y[j]),Sigma)
	PCAQstat = min(z.WI,z.Wald,z.VC)
    z[i,j] = mean(PCAQ_ref < PCAQstat)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PCAQ",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()


#############
C = 0.05
n=100
x = seq(-5,5,length=n)
y = seq(-5,5,length=n)
J = as.matrix(c(1,1))
a = t(J)%*%solve(Sigma)
z = matrix(NA,ncol=n,nrow=n)
PCMeta_ref = PCMetaEmpDist(simNum=1e5, Sigma=Sigma)
for(i in 1:n){
   for(j in 1:n){
    z.PCMinP = PCMinP(c(x[i],y[j]),Sigma)
	z.PCFisher = PCFisher(c(x[i],y[j]),Sigma)
	z.PCLC = PCLC(c(x[i],y[j]),Sigma)
    z.WI = WI(c(x[i],y[j]),Sigma)
	z.Wald = Wald(c(x[i],y[j]),Sigma)
	z.VC= VC(c(x[i],y[j]),Sigma)
	PCMetastat = min(z.PCMinP,z.PCFisher,z.PCLC, z.WI,z.Wald,z.VC)
    z[i,j] = mean(PCMeta_ref < PCMetastat)
   }
}
contour(x,y,z,levels=C,drawlabels=FALSE,col="black",font.lab=font.lab,cex.axis=cex.axis,lty=lty,cex.lab=cex.lab,lwd=lwd,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
title( "PCO",cex.main=cex.main,font.main=font.main)
addPoint()
addLines()
addArrows()

# C = 0.05
# n=100
# x = seq(-5,5,length=n)
# y = seq(-5,5,length=n)
# J = as.matrix(c(1,1))
# a = t(J)%*%solve(Sigma)
# z = matrix(NA,ncol=n,nrow=n)
# for(i in 1:n){
   # for(j in 1:n){
    # z[i,j] = MinP(c(x[i],y[j]),Sigma)
   # }
# }
# contour(x,y,z,levels=C,drawlabels=FALSE,col="black",cex.lab=0.85,xlim=c(-6,6),ylim=c(-6,6),xlab ="Z1",ylab="Z2" )
# title( "MinP",cex=2)
# addPoint()
# addLines()
# addArrows()

dev.off()
