library(mvtnorm)
library(CompQuadForm)
library(MPAT)
source("./tatesLiu.R")
# Accept UNIX command line arguments
args=(commandArgs(TRUE))
if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

####### input file 

datapath = "~/zliulab/Rany/ResultPCAJASAGC/"
inputfile = paste0(datapath,"M10Zscores.txt")
dat = read.table(inputfile,header=T)

## output file 
outDir = "~/zliulab/Rany/ResultPCAJASAGC/"
modName = paste0("M",as.character(mod1))
outputFileName = paste0(outDir,modName,"ResultnoPCOUnitvar",".txt.gz")
outputFileName

### take the Z-score columns only
### based on qq plot of first round analysis 
gc.BMI = 1.098555
gc.DBP = 1.024134
gc.fGLU = 1.067385
gc.fINSULIN = 1.073763
gc.HDL = 1.011243
gc.LDL = 1.010773
gc.SBP = 1.019964
gc.TC = 1.0009363
gc.TG = 1.004205
gc.WHR = 1.047457

### gc correction by sqrt(lambda)
dat$Z.BMI = dat$Z.BMI/sqrt(gc.BMI)
dat$Z.DBP = dat$Z.DBP/sqrt(gc.DBP)
dat$Z.fGLU = dat$Z.fGLU/sqrt(gc.fGLU)
dat$Z.fINSULIN = dat$Z.fINSULIN/sqrt(gc.fINSULIN)
dat$Z.HDL = dat$Z.HDL/sqrt(gc.HDL)
dat$Z.LDL = dat$Z.LDL /sqrt(gc.LDL)
dat$Z.SBP = dat$Z.SBP/sqrt(gc.SBP)
dat$Z.TC = dat$Z.TC/sqrt(gc.TC)
dat$Z.TG = dat$Z.TG/sqrt(gc.TG)
dat$Z.WHR = dat$Z.WHR/sqrt(gc.WHR)

############################
ZCols = grep("Z.",names(dat),value=T) 
ZCols
### if mod1==8, remove TC, DBP
if(mod1==8){
    ZCols = ZCols[-which(ZCols=="Z.TC"|ZCols=="Z.DBP")]
}
####
ZCols 

Z.mat = dat[,ZCols]

summary(Z.mat)
head(Z.mat)
dim(Z.mat)
snp.num = dim(Z.mat)[1]
snp.num
## must transform to a matrix!!!!!!!! 
Z.mat = as.matrix(Z.mat[complete.cases(Z.mat), ]) ## remove NA values again!! 
snp.num = dim(Z.mat)[1]
snp.num
head(Z.mat)
dim(Z.mat)
cov(Z.mat)
Z.mat = scale(Z.mat,center=F,scale=T)


######################################################
## use all the SNPs, results similar to LD pruning 
#####################################################
Sigma = cov(Z.mat)
Sigma
cor(Z.mat)
#######

p.PC1 = rep(NA,snp.num)
p.PC2 = rep(NA,snp.num)
p.PC3 = rep(NA,snp.num)
p.PC4 = rep(NA,snp.num)
p.PC5 = rep(NA,snp.num)
p.PC6 = rep(NA,snp.num)
p.PC7 = rep(NA,snp.num)
p.PC8 = rep(NA,snp.num)
#p.PC9 = rep(NA,snp.num)
#p.PC10 = rep(NA,snp.num)


p.PCLC = rep(NA,snp.num)
p.PCFisher = rep(NA,snp.num)
p.PCMinP = rep(NA,snp.num)

p.WI = rep(NA,snp.num)
p.Wald = rep(NA,snp.num)
p.VC = rep(NA,snp.num)
p.tates = rep(NA,snp.num)


### number of phenotypes 
nvar = dim(Sigma)[1]
nvar 
### this code can only handle nvar=9 or 10 
## loop over all snps
for(i in 1:snp.num){
## loop over all PCs
Z = Z.mat[i,]
p.PC1[i] = PC(Z,Sigma,1)
p.PC2[i] = PC(Z,Sigma,2)
p.PC3[i] = PC(Z,Sigma,3)
p.PC4[i] = PC(Z,Sigma,4)
p.PC5[i] = PC(Z,Sigma,5)
p.PC6[i] = PC(Z,Sigma,6)
p.PC7[i] = PC(Z,Sigma,7)
p.PC8[i] = PC(Z,Sigma,8)
#p.PC9[i] = PC(Z,Sigma,9)
#if(nvar==10){
#p.PC10[i] = PC(Z,Sigma,10)
#}

p.PCLC[i] = PCLC(Z,Sigma)
p.PCFisher[i] = PCFisher(Z,Sigma)
p.PCMinP[i] = PCMinP(Z,Sigma)
p.WI[i] = WI(Z,Sigma)
p.Wald[i] = Wald(Z,Sigma)
p.VC[i] = VC(Z,Sigma)

p.tates[i] = tatesLiu(Z,Sigma)
if(i%%10000==0){
print(i)}
}

#######################################
# output to a gz file to save space 
#######################################
dat$p.PC1 = p.PC1
dat$p.PC2 = p.PC2
dat$p.PC3 = p.PC3
dat$p.PC4 = p.PC4
dat$p.PC5 = p.PC5
dat$p.PC6 = p.PC6
dat$p.PC7 = p.PC7
dat$p.PC8 = p.PC8
#dat$p.PC9 = p.PC9
#if(nvar==10){
#dat$p.PC10 = p.PC10
#}

dat$p.PCLC = p.PCLC
dat$p.PCFisher = p.PCFisher
dat$p.PCMinP = p.PCMinP

dat$p.WI = p.WI
dat$p.Wald = p.Wald
dat$p.VC = p.VC
dat$p.tates = p.tates


head(dat)
dim(dat)


########
write.table(dat,gzfile(outputFileName),sep="\t",quote=F,row.names=F)



