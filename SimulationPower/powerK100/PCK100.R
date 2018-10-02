
## PCMinP simulation and compare powers with WI, Wald and VC
rm(list=ls())
ls()
library(mvtnorm)
library(MPAT)
#library(GHC)
source("./tatesLiu.R")
## principal angles
PAfun = function(Z,u){
  a = sum(Z*u)/(sqrt(sum(Z^2)*sum(u^2)))
  b=acos(a)*180/pi
  b 
}

l2norm = function(a.vec){
sqrt(sum(a.vec^2))
}

Oracle = function(mu,Z.vec,sigma){
  oracle.test = (t(mu)%*%solve(sigma)%*%Z.vec)^2/((t(mu)%*%solve(sigma)%*%mu))
  pchisq(oracle.test,df=1,lower.tail=F)
}

##############33
##################################
PCAQEmpDist = function(simNum=1e5,sigma=sigma){
  # Step 1: simulation Z-vectors under the null
  # equivalent to parametric Bootstrap
  K = dim(sigma)[1]
  ## this multivariate normal must be under the Null!!
  ## This is the basis for generating a null reference distribution
  Z.mat = rmvnorm(simNum,rep(0,K),sigma = sigma)
  PCQMinP = rep(NA,simNum)
  for(b in 1:simNum){
    p.WI = WI(Z.mat[b,],sigma)
    p.Wald = Wald(Z.mat[b,],sigma)
    p.VC = VC(Z.mat[b,],sigma)
    ## compute WI, Wald, VC, take minimum
    PCQMinP[b] = min(c(p.WI,p.Wald,p.VC))
  }
  ## return the empirical null distribution
  return(PCQMinP)
}


PCMetaEmpDist = function(simNum=1e5,sigma=sigma){
  # Step 1: simulation Z-vectors under the null
  # equivalent to parametric Bootstrap
  K = dim(sigma)[1]
  ## this multivariate normal must be under the Null!!
  ## This is the basis for generating a null reference distribution
  Z.mat = rmvnorm(simNum,rep(0,K),sigma = sigma)
  PCMetaMinP = rep(NA,simNum)
  for(b in 1:simNum){
    p.PCMinP = PCMinP(Z.mat[b,],sigma)
    p.PCFisher = PCFisher(Z.mat[b,],sigma)
    p.PCLC = PCLC(Z.mat[b,],sigma)
    p.WI = WI(Z.mat[b,],sigma)
    p.Wald = Wald(Z.mat[b,],sigma)
    p.VC = VC(Z.mat[b,],sigma)
    ## compute WI, Wald, VC, take minimum
    PCMetaMinP[b] = min(c(p.PCMinP,p.PCFisher,p.PCLC,p.WI,p.Wald,p.VC))
  }
  ## return the empirical null distribution
  return(PCMetaMinP)
}

## Generate a null distribution for PCQ.MinP
## This step is most time-consuming
## we use this null to compute the p-value of PCQMinP
sigma = read.csv("./CorMatK100.csv",header=T)
sigma = as.matrix(sigma)
sigma = round(sigma,4)

t1 = proc.time()
PCAQ.Null = PCAQEmpDist(1e5,sigma = sigma)
print(proc.time() - t1)
print("Done PCAQ Null ")
t2 = proc.time()
PCMeta.Null = PCMetaEmpDist(1e5,sigma = sigma)
print(proc.time()-t2)
print("Done PCO Null ")
#######################################
### simulation function 
############################################

PCA.power = function(mu,sigma,sim.n.power = 500,alpha = 0.05,pcoracle=100){
################################################################
## under the alternatives 
###############################################
#############
print("staring")
Z.mat = rmvnorm(sim.n.power,mean=mu,sigma=sigma)
##############################

#############################################
### power comparison under the alternatives 
##############################################
p.Oracle = rep(NA,sim.n.power)
p.pc1 = rep(NA,sim.n.power)
p.pc2 = rep(NA,sim.n.power)
p.pc3 = rep(NA,sim.n.power)
p.MinP = rep(NA,sim.n.power)
p.PCMinP = rep(NA,sim.n.power)
p.PCFisher = rep(NA,sim.n.power)
p.PCLC = rep(NA,sim.n.power)
p.WI = rep(NA,sim.n.power)
p.Wald = rep(NA,sim.n.power)
p.VC = rep(NA,sim.n.power)
p.PCAQ = rep(NA,sim.n.power)
p.PCMeta = rep(NA,sim.n.power)
p.tates = rep(NA,sim.n.power)

print("Looping")
for(i in 1:sim.n.power){
#print(i)
p.Oracle[i] = Oracle(mu=mu,Z.vec = Z.mat[i,],sigma=sigma)
p.pc1[i] = PC(Z.mat[i,],sigma,PCorder=1)
p.pc2[i] = PC(Z.mat[i,],sigma,PCorder=2)
#####
p.pc3[i] = PC(Z.mat[i,],sigma,PCorder=pcoracle)
p.MinP[i] = MinP(Z.mat[i,],sigma)
p.PCMinP[i] = PCMinP(Z.mat[i,],sigma)
p.PCFisher[i] = PCFisher(Z.mat[i,],sigma)
p.PCLC[i] = PCLC(Z.mat[i,],sigma)
p.WI[i] = WI(Z.mat[i,],sigma)
p.Wald[i] = Wald(Z.mat[i,],sigma)
p.VC[i] = VC(Z.mat[i,],sigma)
PCQMinP = min(c(p.WI[i],p.Wald[i],p.VC[i]))
p.PCAQ[i] = mean(PCAQ.Null < PCQMinP)

PCMetaMinP = min(c(p.PCMinP[i],p.PCFisher[i],p.PCLC[i], p.WI[i],p.Wald[i],p.VC[i]))
p.PCMeta[i] = mean(PCMeta.Null < PCMetaMinP )
p.tates[i] = tatesLiu(Z.mat[i,],sigma)
}
### powers at level alpha
pw.Oracle = mean(p.Oracle < alpha)
pw.pc1 = mean(p.pc1 < alpha)
pw.pc2 = mean(p.pc2 < alpha)
pw.pc3 = mean(p.pc3 < alpha)
pw.MinP = mean(p.MinP < alpha)
pw.PCMinP = mean(p.PCMinP < alpha)
pw.PCFisher = mean(p.PCFisher < alpha)
pw.PCLC = mean(p.PCLC < alpha)
pw.WI = mean(p.WI < alpha)
pw.Wald = mean(p.Wald < alpha)
pw.VC = mean(p.VC < alpha)
pw.PCAQ = mean(p.PCAQ < alpha)
pw.PCO = mean(p.PCMeta < alpha)
pw.tates = mean(p.tates <alpha)

output = c(pw.Oracle,pw.pc1,pw.pc2,pw.pc3,pw.MinP,pw.PCMinP,pw.PCFisher,pw.PCLC,pw.WI,pw.Wald,pw.VC,pw.PCAQ,pw.PCO,pw.tates)
return(output) 
}


###########################
## wrap up funciton to run simulation 
########################
runSim = function(mu,sigma,sim.n.power =500,pcoracle=100){
out =  PCA.power(mu,sigma,sim.n.power = sim.n.power,pcoracle = pcoracle) 

U = round(eigen(sigma)$vectors,2)
theta = rep(NA,3) 
K=length(mu)
theta[1] = PAfun(mu,U[,1])
theta[2] = PAfun(mu,U[,2])
theta[3] = PAfun(mu,U[,K])

theta = round(theta,1)
mu.norm = round(l2norm(mu),2)
res = c(mu.norm,theta,out)
res
}



###############################

sigma = read.csv("./CorMatK100.csv",header=T)
sigma = as.matrix(sigma)
sigma = round(sigma,4)
K = dim(sigma)[1]
K

###################################
### spase
eigen(sigma)$values
U = eigen(sigma)$vectors
mu1 = round(U[,K]*0.18,2)
mu1
out1 = runSim(mu1,sigma,sim.n.power=1000,pcoracle = 100)
out1 

## sparse
mu2 = round(U[,95],2)
mu2[11:100]=0
mu2
out2 = runSim(mu2,sigma,sim.n.power=1000,pcoracle=100)
out2

## dense signals
mu3 = round(U[,50],2)*3.4
mu3
out3 = runSim(mu3,sigma,sim.n.power=1000,pcoracle=100)
out3


### PC1, dense        
mu4 = U[,1]*15
mu4
out4 = runSim(mu4,sigma,sim.n.power=1000,pcoracle=100)
out4

## PC80, dense
mu5 = round(U[,80],2)*2
mu5
out5 = runSim(mu5,sigma,sim.n.power=1000,pcoracle=100)
out5

############
result = rbind(out1,out2,out3,out4,out5)
colnames(result)=c("norm","theta1","theta2","thetaK","Oracle","pc1","pc2","pcK","MinP","PCMinP","PCFisher","PCLC","WI","Wald","VC","PCAQ","PCO","TATES")
result 


write.csv(result,file="PCA_powerK100JASA.csv",quote = TRUE) 

betaVec = c(toString(mu1),toString(mu2),toString(mu3),toString(mu4),toString(mu5))
betaVec
result
result.new = cbind(result,betaVec)
result.new
write.csv(result.new,file="PCA_powerK100JASA_meanvectors.csv",quote = TRUE) 
