  ## power simulation funciton 
  ## input: K: number of phenotypes;
  ## input Sigma: correlation structure 
  ## simNum: number of replications 
  ## alpha: significance level 
  ###########################################################
# Estimate the SigmaX for PCAQ p-value computation
# This is a pre-computation step for PCAQ 
# input: the correlation matrix among Z-scores
# output: the SigmaX matrix 
SigmaXEstimate = function(Sigma,simNum=1000){

X.WI = rep(NA,simNum)
X.Wald =rep(NA,simNum)
X.VC =rep(NA,simNum)
for(i in 1:simNum){
K = dim(Sigma)[1]
Z.vec = t(rmvnorm(1,rep(0,K),Sigma))
X.WI[i] = qnorm(WI(Z.vec,Sigma))
X.Wald[i] = qnorm(Wald(Z.vec,Sigma))
X.VC[i] = qnorm(VC(Z.vec,Sigma))
}
## if VC gives pvalue=0, then X.VC is inf, then cor function can't compute it 
# remove all rows with non-finite values, this could happen when X.WI is Inf...
X.mat = cbind(X.WI,X.Wald,X.VC)
X.mat = X.mat[!rowSums(!is.finite(X.mat)),]
SigmaX = cor(X.mat)
return(SigmaX)
}

# Estimate the SigmaMeta for PCMeta p-value computation
# This is a pre-computation step for PCMeta
# input: the correlation matrix among Z-scores
# output: the SigmaMeta matrix 
SigmaMetaEstimate = function(Sigma,simNum=1000){

X.PCMinP = rep(NA,simNum)
X.PCFisher = rep(NA,simNum)
X.PCLC = rep(NA,simNum)
X.WI = rep(NA,simNum)
X.Wald =rep(NA,simNum)
X.VC =rep(NA,simNum)

for(i in 1:simNum){
K = dim(Sigma)[1]
Z.vec = t(rmvnorm(1,rep(0,K),Sigma))

X.PCMinP[i] = qnorm(PCMinP(Z.vec,Sigma))
X.PCFisher[i] = qnorm(PCFisher(Z.vec,Sigma))
X.PCLC[i] = qnorm(PCLC(Z.vec,Sigma))
X.WI[i] = qnorm(WI(Z.vec,Sigma))
X.Wald[i] = qnorm(Wald(Z.vec,Sigma))
X.VC[i] = qnorm(VC(Z.vec,Sigma))
}

X.mat = cbind(X.PCMinP,X.PCFisher,X.PCLC,X.WI,X.Wald,X.VC)
X.mat = X.mat[!rowSums(!is.finite(X.mat)),] ## remove rows containing NA, NaN,Inf
SigmaMeta = cor(X.mat)
return(SigmaMeta)

}
  
  
## input: Z-score: Z.vec; Sigma: the correlation matrix among Z-scores
## SigmaX: is the correlation matrix among X, which is inverse normal transformatoin of MinP 
## SigmaX: should be pre-computed using simulation 
PCAQ = function(Z.vec,Sigma,SigmaX){
 ## avoid Z.vec is a row vector which won't work later
 ## Z.vec should be a column vector object, not a 1*K matrix 
 if(!is.vector(Z.vec)){
 Z.vec = as.vector(Z.vec)
 }


p.WI = WI(Z.vec,Sigma)
p.Wald = Wald(Z.vec,Sigma)
p.VC = VC(Z.vec,Sigma)
##take minimum
T.minp = min(c(p.WI,p.Wald,p.VC))
## use multivariate normal CDF to compupte p-value

K.X = dim(SigmaX)[1]
if(K.X!=3){
stop("The dimension of SigmaX is not 3!!")
}

lower_bound = rep(qnorm(T.minp), K.X)
upper_bound = rep(Inf,K.X)
a = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
        K.X), sigma = SigmaX, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
p.PCAQ = 1 - a[1]
return(p.PCAQ)
}



## input: Z-score: Z.vec; Sigma: the correlation matrix among Z-scores
## SigmaMeta: is the correlation matrix among X, which is inverse normal transformatoin of MinP 
## SigmaMeta: should be pre-computed using simulation 
PCMeta = function(Z.vec,Sigma,SigmaMeta){
 ## avoid Z.vec is a row vector which won't work later
 ## Z.vec should be a column vector object, not a 1*K matrix 
 if(!is.vector(Z.vec)){
 Z.vec = as.vector(Z.vec)
 }

##########
p.PCMinP = PCMinP(Z.vec,Sigma)
p.PCFisher = PCFisher(Z.vec,Sigma)
p.PCLC = PCLC(Z.vec,Sigma)
p.WI = WI(Z.vec,Sigma)
p.Wald = Wald(Z.vec,Sigma)
p.VC = VC(Z.vec,Sigma)
T.minp = min(p.PCMinP,p.PCFisher,p.PCLC, p.WI,p.Wald,p.VC)

K.X = dim(SigmaMeta)[1]
if(K.X!=6){
stop("The dimension of SigmaMeta is not 6!!")
}
lower_bound = rep(qnorm(T.minp), K.X)
upper_bound = rep(Inf,K.X)

a = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
        K.X), sigma = SigmaMeta, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
pval = 1 - a[1]
return(pval)
}


PowerSim = function(Sigma = Sigma, simNum = 1e4,alpha){
  K = dim(Sigma)[1]
  pval.PCMinP = rep(NA,simNum)
  pval.PCFisher = rep(NA,simNum)
  pval.PCLC = rep(NA,simNum)
  pval.WI = rep(NA,simNum)
  pval.VC = rep(NA,simNum)
  pval.Wald = rep(NA,simNum)
  pval.PCAQ = rep(NA,simNum)
  pval.PCMeta = rep(NA,simNum)
  SigmaX = SigmaXEstimate(Sigma,1e4)
  SigmaMeta = SigmaMetaEstimate(Sigma,1e4)
   
  for(i in 1:simNum){
  Z.vec = rmvnorm(n=1,rep(0,K),sigma = Sigma)
  Z.vec = t(Z.vec) ## want a column vector
  pval.PCMinP[i] = PCMinP(Z.vec,Sigma=Sigma)
  pval.PCFisher[i]=PCFisher(Z.vec,Sigma=Sigma)
  pval.PCLC[i] = PCLC(Z.vec,Sigma=Sigma)
  pval.WI[i] = WI(Z.vec,Sigma=Sigma)
  pval.Wald[i] = Wald(Z.vec, Sigma = Sigma)
  pval.VC[i] = VC(Z.vec,Sigma=Sigma)
  
    Tminp.PCAQ = min(c(pval.WI[i],pval.Wald[i],pval.VC[i])) 

	
	lower_bound = rep(qnorm(Tminp.PCAQ), 3)
	upper_bound = rep(Inf,3)
	a = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
			3), sigma = SigmaX, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
	pval.PCAQ[i] = 1 - a[1]
	
	Tminp.PCMeta = min(c(pval.PCMinP[i],pval.PCFisher[i],pval.PCLC[i],pval.WI[i],pval.Wald[i],pval.VC[i]))	
	upper_bound = rep(Inf,6)
	lower_bound = rep(qnorm(Tminp.PCMeta), 6)

	b = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
			6), sigma = SigmaMeta, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
	pval.PCMeta[i] = 1 - b[1]
	  
  }
  ## count the number of rejections
  power.PCMinP = sum(pval.PCMinP < alpha)
  power.PCFisher = sum(pval.PCFisher < alpha)
  power.PCLC = sum(pval.PCLC < alpha)
  power.WI = sum(pval.WI < alpha)
  power.Wald = sum(pval.Wald < alpha)
  power.VC = sum(pval.VC<alpha)
  power.PCAQ = sum(pval.PCAQ < alpha)
  power.PCMeta = sum(pval.PCMeta < alpha)
  
  res = c(power.PCMinP,power.PCFisher,power.PCLC,power.WI,power.Wald,power.VC,power.PCAQ,power.PCMeta)
  #names(res)=c("PCMinP","PCFisher","PCLC","WI","Wald","VC","PCAQ","PCMeta")
  return(res)
}