##############################
## Rotation plot
###############################
dat = read.csv("PC_Rotation12tests_rho06.csv",header=T)
head(dat)

### test is a vector of test names 
plotRotation = function(powerMat,test,titleName,rho=0.6,color,line.mark=c(135,225)){
p = dim(powerMat)[2]
thetaV = powerMat[,1]
Powers = powerMat[,test]
matplot(thetaV,Powers,type="l", lwd=2,lty=1:dim(Powers)[2],col=color,ylim=c(0,1),xaxt="n",xlab=bquote(paste("Rotation Angle ",phi)), ylab="Power",cex.lab=1.5)
#legend("topleft", legend = test, lty=1:p,col=1:p,cex=0.8) 
legend("topleft", legend = test,col=color,lty=1:dim(Powers)[2],cex=1.5)
abline(v=line.mark ,lty=3,cex=1.5,col="blue")
x.tick = c(0,45,90,135,180,225,270,315,360)
axis(1, at=x.tick, labels=as.character(x.tick))
title(titleName,cex.main=1.8,font.main=1)
}

##############################################
## plotting
############################################


pdf("RotationV2Blackwhite.pdf",height=8,width=18)
par(mfrow=c(2,3))
plotRotation(powerMat = dat,test=c("Oracle","PC1","PC2"),titleName="PC1 versus PC2",color=c(2,3,4))
plotRotation(powerMat = dat,test=c("Oracle","PC1","WI"),titleName="PC1 versus WI",color=c(2,3,5))
plotRotation(powerMat = dat,test=c("Oracle","PC2","VC"),titleName="PC2 versus VC",color=c(2,4,"black"))
plotRotation(powerMat = dat,test=c("Oracle","WI","Wald","VC"),titleName="PCQ",color=c(2,5,"blueviolet","black"))
plotRotation(powerMat = dat,test=c("Oracle","PCFisher","Wald","PCLC"),titleName="Linear versus Nonlinear",color=c(2,"chartreuse2","blueviolet","gray8"),line.mark=c(90,270))
plotRotation(powerMat = dat,test=c("Oracle","PCMinP","PCAQ","PCO"),titleName="Adaptive",color=c(2,"magenta1","slateblue4","cyan2"))
dev.off()