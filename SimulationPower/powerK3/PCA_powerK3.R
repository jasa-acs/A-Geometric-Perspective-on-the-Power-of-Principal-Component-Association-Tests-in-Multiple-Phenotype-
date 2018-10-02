## PCMinP simulation and compare powers with WI, Wald and VC
rm(list=ls())
ls()
library(mvtnorm)
library(MPAT)
#install.packages("pracma")
library(pracma) ## for gram schmidt
source("./tatesLiu.R")


## principal angles
PAfun = function(Z,u){
  Z = Z/sqrt(sum(Z^2))
  u = u/sqrt(sum(u^2))
  a = sum(Z*u)
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


#######################################
### simulation function 
############################################

PCA.power = function(mu,sigma,sim.n.power = 500,alpha = 0.05){
################################################################
## under the alternatives 
###############################################
#############
Z.mat = rmvnorm(sim.n.power,mean=mu,sigma=sigma)
##############################


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
PCAQ.Null = PCAQEmpDist(1e5,sigma = sigma)
PCMeta.Null = PCMetaEmpDist(1e5,sigma = sigma)

#############################################
### power comparison under the alternatives 
##############################################
p.Oracle = rep(NA,sim.n.power)
p.pc1 = rep(NA,sim.n.power)
p.pc2 = rep(NA,sim.n.power)
p.pc3 = rep(NA,sim.n.power)
#p.pc4 = rep(NA,sim.n.power)
p.PCMinP = rep(NA,sim.n.power)
p.PCFisher = rep(NA,sim.n.power)
p.PCLC = rep(NA,sim.n.power)
p.WI = rep(NA,sim.n.power)
p.Wald = rep(NA,sim.n.power)
p.VC = rep(NA,sim.n.power)
p.PCAQ = rep(NA,sim.n.power)
p.PCMeta = rep(NA,sim.n.power)
p.tates = rep(NA,sim.n.power)

for(i in 1:sim.n.power){
p.Oracle[i] = Oracle(mu=mu,Z.vec = Z.mat[i,],sigma=sigma)
p.pc1[i] = PC(Z.mat[i,],sigma,PCorder=1)
p.pc2[i] = PC(Z.mat[i,],sigma,PCorder=2)
p.pc3[i] = PC(Z.mat[i,],sigma,PCorder=3)
#p.pc4[i] = PC(Z.mat[i,],sigma,PCorder=4)
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
#pw.pc4 = mean(p.pc4 < alpha)
pw.PCMinP = mean(p.PCMinP < alpha)
pw.PCFisher = mean(p.PCFisher < alpha)
pw.PCLC = mean(p.PCLC < alpha)
pw.WI = mean(p.WI < alpha)
pw.Wald = mean(p.Wald < alpha)
pw.VC = mean(p.VC < alpha)
pw.PCAQ = mean(p.PCAQ < alpha)
pw.PCMeta = mean(p.PCMeta < alpha)
pw.tates = mean(p.tates <alpha)


output = c(pw.Oracle,pw.pc1,pw.pc2,pw.pc3,pw.PCMinP,pw.PCFisher,pw.PCLC,pw.WI,pw.Wald,pw.VC,pw.PCAQ,pw.PCMeta,pw.tates)
return(output) 
}

###################################
# unstructured correlation 
######################################
K=3
sigma.vec = c(1.00000000, -0.07663404, 0.1553153, -0.4222786,
        -0.07663404,  1.00000000, 0.8787142,  0.2708553,
        0.15531530,  0.87871424, 1.0000000,  0.3842118,
        -0.42227857,  0.27085532, 0.3842118,  1.0000000)
sigma2 = matrix(round(sigma.vec,2),ncol=4,byrow=T)
sigma2 = sigma2[-2,-2]

U = round(eigen(sigma2)$vectors,2)
lambdas = round(eigen(sigma2)$values,2)
lambdas
sum(lambdas)

###########################
## wrap up funciton to run simulation 
########################
runSim = function(mu,sigma2,sim.n.power =500){
#U = round(eigen(sigma2)$vectors,2)
U = eigen(sigma2)$vectors
out =  PCA.power(mu,sigma2,sim.n.power = sim.n.power) 
theta = rep(NA,K)
for(i in 1:3){
  theta[i] = PAfun(mu,U[,i])
}
mu = round(mu,2)
theta = round(theta,1)
mu.norm = round(l2norm(mu),1)
res = c(toString(mu),mu.norm,theta,out)
res
}

#################
## 1st PC 
################
mu1= U[,1]*3.8
out1 = runSim(mu1,sigma2)
out1
#################
## 2nd PC 
################
mu2 = U[,2]*3.5
out2 = runSim(mu2,sigma2)
out2 

##############
## 3nd PC 
#############
mu3 = U[,3]*1.8
mu3
out3 = runSim(mu3,sigma2)
out3 


#####################################
### in the middle of four PCs
### mu5 has equal angles with all four PCs
### PCLC is Oracle 
###############################
mu5 = round(rowMeans(U)*4,2)
mu5
round(U%*% rep(1,K),2) 
### check the angles 
PAfun(mu5,U[,1])
PAfun(mu5,U[,2])
PAfun(mu5,U[,3])
## run simulation 
out5 = runSim(mu5,sigma2)
out5



############################################################
### mu6 PCLC acheives its maximal power, not Oracle 
###########################################################
U = eigen(sigma2)$vectors 
 lambdas = eigen(sigma2)$values
 coef.beta = (1/lambdas)/sqrt(sum(1/lambdas^2))
 mu6  = U %*%coef.beta*0.7 
 mu6
# #######
# Lambda.mat = diag(lambdas,nrow=K)
# beta.direction = U %*%solve(Lambda.mat)%*%rep(1,K)
# beta.direction/mu6
# ##########
# ### check the angles 
# PAfun(mu6,U[,1])
# PAfun(mu6,U[,2])
# PAfun(mu6,U[,3])
# ## run simulation 
# out6 = runSim(mu6,sigma2)
# out6 

### mu8 PCLC least powerful, powerless.  
#sin.k = (1/lambdas)/sqrt(sum(1/lambdas^2))
#coef.beta = sqrt(1-sin.k^2) ## positive or negative cos function?? never  know, this is wrong 
A =cbind(mu6,rep(1,K)) 
Q = gramSchmidt(A)$Q
Q[,1]/mu6 ## proportional 
mu7 = round(Q[,2],2)*3.3
sum(mu6*mu7) ## check orthogonal 
mu7
out7 = runSim(mu7,sigma2)
out7

##################################
## equal powers for all PCs
#################################
# mu6 equal powers for all PCs, PCMinP will be least powerful?? not really 
# U = eigen(sigma2)$vectors 
# coef.beta = sqrt(lambdas/K)
# coef.beta
# mu8 =  U %*%coef.beta*3.1
# mu8
# ### check the angles 
# theta1 = PAfun(mu8,U[,1])
# theta2 = PAfun(mu8,U[,2])
# theta3 = PAfun(mu8,U[,3])
# 
# theta1
# theta2
# theta3
# 
# ## check ncp for the four PCs, these 4 ncps should be equal 
# (cos(theta1*pi/180))^2/lambdas[1]
# (cos(theta2*pi/180))^2/lambdas[2]
# (cos(theta3*pi/180))^2/lambdas[3]
# 
# out8 = runSim(mu8,sigma2)
# out8

##############
### no PC4, so no out4
result = rbind(out1,out2,out3,out5,out7)
colnames(result)=c("mean","norm","theta1","theta2","thetaK","Oracle","pc1","pc2","pcK","PCMinP","PCFisher","PCLC","WI","Wald","VC","PCAQ","PCO","TATES")
result
## the mean is correctly outputed  when quote is true
write.csv(result,file="PCA_powerK3Jasa.csv",quote = TRUE) 





