
# Accept UNIX command line arguments
args=(commandArgs(trailingOnly=TRUE))
if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}


library(mvtnorm)
library(MPAT)
source("PowerSim.R")
##receive unix command line arguments

simNum = 100
alpha =0.05
ncore = 4

#################################
## unstr means: Sigma is unstructured or not 

Sigma = read.csv("CorMatK100.csv",header=T)
Sigma = as.matrix(Sigma)
Sigma 
K = dim(Sigma)[1]
###############################


###############################
## multiple core computing
###############################
library(doMC)
ncore = ncore
ncore = 4
registerDoMC(ncore)


ptm = proc.time()
res = foreach(icount(ncore),.combine='rbind') %dopar%{
res = PowerSim(Sigma=Sigma,simNum = simNum,alpha=alpha)
res
}
proc.time() - ptm 
res

ncore = dim(res)[1]
ncore
res = as.data.frame( t(colSums(res)/(ncore*simNum)))
res

resdir="./Result/"
filename = paste0(resdir,"K_",as.character(K),"alpha",as.character(alpha),".csv")


names(res) = c("PCMinP","PCFisher","PCLC","WI","Wald","VC","PCAQ","PCO")
res
write.csv(res,file=filename,quote=F,row.names=F)
