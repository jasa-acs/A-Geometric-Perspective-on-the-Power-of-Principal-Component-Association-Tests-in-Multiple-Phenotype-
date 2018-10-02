###################################
# Author: Zhonghua Liu
# Date: 09/24/2015
# For Metbolic syndrome with Rany
###################################
library(mvtnorm)
library(CompQuadForm)
library(MPAT)
source("./tatesLiu.R")
SBP_in = FALSE
DBP_in = FALSE
BMI_in = FALSE
WHR_in = FALSE
fGLU_in = FALSE
fINSULIN_in = FALSE
LDL_in = FALSE
TG_in = FALSE
HDL_in = FALSE
TC_in = FALSE
# Accept UNIX command line arguments
args=(commandArgs(TRUE))
if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}
## output file 
outDir = "~/zliulab/Rany/ResultPCAJASA/"
modName = paste0("M",as.character(mod1),"_",as.character(mod2))
outputFileName = paste0(outDir,modName,".txt.gz")
outputFileName

#source("pickSNPs.R") ## used only for LD pruning..almost the same as using all the SNPs
source("DataPrep.R")
source("DataLoad.R")
source("multmerge.R")


#################################################
## merge files based on the chr:pos column 
###################################################
##check the objects in current environment 
ls()
## return a list of objects with names "dat.*"
datalist = mget(grep("dat.",ls(),value=T))
datalist
dat.merge = multmerge(datalist)
head(dat.merge)
### take the Z-score columns only 
Z.mat = dat.merge[,grep("Z.",names(dat.merge))]
summary(Z.mat)
snp.num = dim(Z.mat)[1]
snp.num
## must transform to a matrix!!!!!!!! 
Z.mat = as.matrix(Z.mat[complete.cases(Z.mat), ]) ## remove NA values again!! 
snp.num = dim(Z.mat)[1]
snp.num
head(Z.mat)
######################################################
## use all the SNPs, results similar to LD pruning 
#####################################################
Sigma = cor(Z.mat)
Sigma
write.csv(Sigma,file=paste0("./log/",modName,"Sigma.csv"))
source("DataAnalysis.R")
head(dat.merge)
write.table(dat.merge,gzfile(outputFileName),quote=F,row.names=F)



