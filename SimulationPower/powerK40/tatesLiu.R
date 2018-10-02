## this function is to calculate TATES p-value 
## input: Z-score vector, correlation matrix among Z-scores
## output: an overall p-value 
tatesLiu = function(Z.vec,cormat){

pvalues = pchisq(Z.vec^2,df=1,lower.tail=F)
nvar = length(pvalues)
# read in the beta weights from the 6th order polynomial
betaa=c(-0.0007609278,-0.0022965148,0.6226249243,0.0148755138,
         0.1095155903,-0.0218930325,0.2178970393)

ps = pvalues 
tmp=sort(ps,index.return=T)  # gives location of the sorted p-values
pj=tmp$x # sorted p-values
iorder=tmp$ix # index
# 2.----------------- get the correlation matrix between the variables 

r2=matrix(0,nvar,nvar)
# order correlation matrix according to the rankorder of the p-values
## this place use the correlation matrix 
r = cormat 
r2=r[iorder,iorder]

# 3. ---------------- weight symptom correlations [r2] with regression weights 
# matrix ro contains the predicted correlations between p-values,
# i.e., predicted from the correlations between variables
ro=diag(nvar)
for (i1 in 1:nvar)  {  
for (i2 in 1:i1) {
if (i1>i2) {
er=r2[i1,i2]
ro[i1,i2]=ro[i2,i1]= betaa[7]*er^6+betaa[6]*er^5+betaa[5]*er^4+betaa[4]*er^3+betaa[3]*er^2+betaa[2]*er+betaa[1]
}}}

# 4. ---------------- determine eigen values based on entire p-value matrix 
# get Mall (Me based on all p-values)
alllam=eigen(ro[1:nvar,1:nvar])$values #eigenvalues of the ro matrix
mepj=nvar
for (i1 in 1:nvar) {
mepj=mepj-(alllam[i1]>1)*(alllam[i1]-1) }

# 5. ---------------- determine eigen values top x p-values 
# sellam is eigenvalues based on varying nr of top SNPs

mej=matrix(c(seq(1,nvar,1)),nvar,1,byrow=T)

for (j in 1:nvar) { 
sellam=eigen(ro[1:j,1:j])$values #eigenvalues of the ro matrix
id=j
# subtract 1-eigenvalue for those eigenvalues >1    # page 284 Li et al., 2011
for (i1 in 1:id) {
mej[j,1]=mej[j,1]-(sellam[i1]>1)*(sellam[i1]-1)
}
}

# 6. ---------------- weight sorted p-values with eigenvalues ratio 

pg=matrix(0,nvar,1)
for (i in 1:nvar) {
pg[i,1]=(mepj/mej[i,1])*pj[i]
}

pg=pg[iorder]  # p-values back in original order so that pval[1] corresponds to item[1] etc !!

return(min(pg))
}