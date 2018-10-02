## All the input is a Z column vector with length K 
library(mvtnorm)
library(MPAT)

# Accept UNIX command line arguments
args=(commandArgs(trailingOnly=TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}


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
#####################
########################
## data input 
###########################
## command line argument 
chunkid
datapath = "~/zliulab/Rany/ResultPCAJASAGC/chunkforPCOinput/"
modName = paste0("M",as.character(mod1))
inputfile = paste0(datapath,modName,"chunk",as.character(chunkid),".txt.gz")
dat = read.table(inputfile,header=T)
head(dat)
dim(dat)
### take the Z-score columns only 
snp.num = dim(dat)[1]
snp.num
######################################################
## use all the SNPs, results similar to LD pruning 
#####################################################
Sigma = read.csv(paste0(modName,"Cormat.csv"),header=T)
Sigma = as.matrix(Sigma)
Sigma



SigmaX = SigmaXEstimate(Sigma,1e4)
SigmaX
SigmaMeta = SigmaMetaEstimate(Sigma,1e4)
SigmaMeta 

p.PCAQ = rep(NA,snp.num)
p.PCO = rep(NA,snp.num)

for(i in 1:snp.num){

p.PCMinP = dat$p.PCMinP[i]
p.PCFisher = dat$p.PCFisher[i]
p.PCLC = dat$p.PCLC[i] 
p.WI = dat$p.WI[i]
p.Wald = dat$p.Wald[i]
p.VC = dat$p.VC[i]

#######
Tminp.PCAQ = min(c(p.WI,p.Wald,p.VC)) 
lower_bound = rep(qnorm(Tminp.PCAQ), 3)
upper_bound = rep(Inf,3)
a = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
			3), sigma = SigmaX, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
p.PCAQ[i] = 1 - a[1]
#######	
Tminp.PCMeta = min(c(p.PCMinP,p.PCFisher,p.PCLC,p.WI,p.Wald,p.VC))	
upper_bound = rep(Inf,6)
lower_bound = rep(qnorm(Tminp.PCMeta), 6)
b = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
			6), sigma = SigmaMeta, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
p.PCO[i] = 1 - b[1]

}


####################################
## output file 
####################################
outDir =  "~/zliulab/Rany/ResultPCAJASAGC/chunkResultPCO"
modName = paste0("M",as.character(mod1))
outputFileName = paste0(outDir,modName,"chunkResultPCO",as.character(chunkid),".txt.gz")
outputFileName

dat$p.PCAQ = p.PCAQ
dat$p.PCO = p.PCO



write.table(dat,file=outputFileName,sep="\t",quote=F,row.names=F)













