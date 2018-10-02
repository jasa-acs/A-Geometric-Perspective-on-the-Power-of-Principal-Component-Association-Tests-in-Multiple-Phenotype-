
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
K =3
rho = 0.3
simNum = 100
alpha =0.05
unstr=TRUE
ncore = 4


#################################
## unstr means: Sigma is unstructured or not 
if(unstr == FALSE){
J = rep(1,K)
Sigma = rho * J%*%t(J) + (1-rho)*diag(K)
}

if(unstr==TRUE){
Sigma.vec = c(1.00000000, -0.07663404, 0.1553153, -0.4222786,
        -0.07663404,  1.00000000, 0.8787142,  0.2708553,
        0.15531530,  0.87871424, 1.0000000,  0.3842118,
        -0.42227857,  0.27085532, 0.3842118,  1.0000000)
Sigma = matrix(round(Sigma.vec,2),ncol=4,byrow=T)
Sigma = round(Sigma,2)
}


Sigma = Sigma[-2,-2]
Sigma
###############################

###############################
## multiple core computing
###############################
library(doMC)
ncore = ncore
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


resdir="./Result/"
if(unstr==FALSE){
filename = paste0(resdir,"K_",as.character(K),"_rho_",as.character(rho),"alpha",as.character(alpha),".csv")
}
if(unstr==TRUE){
filename = paste0(resdir,"K_",as.character(K),"_unstr","alpha",as.character(alpha),".csv")
}

names(res) = c("PCMinP","PCFisher","PCLC","WI","Wald","VC","PCAQ","PCO")
write.csv(res,file=filename,quote=F,row.names=F)
