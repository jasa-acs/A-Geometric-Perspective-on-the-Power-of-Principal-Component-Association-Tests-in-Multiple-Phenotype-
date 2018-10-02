p.PCFisher = rep(NA,snp.num)
p.PC1 = rep(NA,snp.num)
p.PC2 = rep(NA,snp.num)
p.PC3 = rep(NA,snp.num)
p.PC4 = rep(NA,snp.num)
p.PC5 = rep(NA,snp.num)
p.PC6 = rep(NA,snp.num)
p.PC7 = rep(NA,snp.num)
p.PC8 = rep(NA,snp.num)
p.PC9 = rep(NA,snp.num)
p.PC10 = rep(NA,snp.num)
p.PCLC = rep(NA,snp.num)
p.WI = rep(NA,snp.num)
p.Wald = rep(NA,snp.num)
p.VC = rep(NA,snp.num)
p.PCMinP = rep(NA,snp.num)
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
p.PC9[i] = PC(Z,Sigma,9)
if(nvar==10){
p.PC10[i] = PC(Z,Sigma,10)
}
p.PCLC[i] = PCLC(Z,Sigma)
p.PCFisher[i] = PCFisher(Z,Sigma)
p.WI[i] = WI(Z,Sigma)
p.Wald[i] = Wald(Z,Sigma)
p.VC[i] = VC(Z,Sigma)
p.PCMinP[i] = PCMinP(Z,Sigma)
p.tates[i] = tatesLiu(Z,Sigma)
}

#######################################
# output to a gz file to save space 
#######################################
dat.merge$p.PC1 = p.PC1
dat.merge$p.PC2 = p.PC2
dat.merge$p.PC3 = p.PC3
dat.merge$p.PC4 = p.PC4
dat.merge$p.PC5 = p.PC5
dat.merge$p.PC6 = p.PC6
dat.merge$p.PC7 = p.PC7
dat.merge$p.PC8 = p.PC8
dat.merge$p.PC9 = p.PC9
if(nvar==10){
dat.merge$p.PC10 = p.PC10
}
dat.merge$p.PCLC = p.PCLC
dat.merge$p.WI = p.WI
dat.merge$p.Wald = p.Wald
dat.merge$p.VC = p.VC
dat.merge$p.PCFisher = p.PCFisher
dat.merge$p.PCMinP = p.PCMinP
dat.merge$p.tates = p.tates
