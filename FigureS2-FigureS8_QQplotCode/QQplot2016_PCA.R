library(GenABEL)
## Accept UNIX command line arguments
args=(commandArgs(trailingOnly=TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

##########################################
### QQ plot function, can't handle p=0 !! 
##########################################
ggd.qqplot = function(pvector, main=NULL, ...) {
    pvector[pvector==0 | pvector < 0 ] = min(pvector[pvector > 0])
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
}

###########################
### load in the results 
#########################
input.path = "~/zliulab/Rany/ResultPCAJASAGC/chunkResultPCOUnitvarM8/M8/"
FileName = paste0("M",as.character(mod1),"chunkResultMergePCOUnitvar.txt")
inputFileName = paste0(input.path ,FileName)
inputFileName
dat = read.table(inputFileName,header=T)
dim(dat)
head(dat)
grep("^p.",colnames(dat),value=T)

############################################
### output top hits by Tests and combined 
############################################
outDir = paste0("~/zliulab/Rany/ResultPCAJASAGC/QQUnitvarM8/",paste0("M",as.character(mod1),"/"))

#######################################
## univariate significance SNPs
######################################
sig = 5e-08
ind.uniSig = which(dat$p.BMI < sig | dat$p.DBP < sig | dat$p.fGLU < sig | dat$p.fINSULIN < sig |dat$p.HDL < sig|dat$p.LDL<sig|dat$p.SBP<sig|dat$p.TC < sig| dat$p.TG<sig|dat$p.WHR<sig )
length(ind.uniSig)

## Remove univariate significant SNPs
#dat = dat[-ind.uniSig,]
##################################

myQQ = function(method)
{
col_p_name = paste0("p.",method)
pvals = dat[,col_p_name]
pvals[pvals==0|pvals < 0] = min(pvals[pvals >0])
png(paste0(outDir, "QQ_",method,".png"))
ggd.qqplot(pvals,main=method)
dev.off()
lam.median = estlambda(pvals,method="median")$estimate
lam.reg = estlambda(pvals)$estimate
return(c(lam.median,lam.reg)) ## lambda values 
}

pvalCols = grep("^p.",colnames(dat),value=T)
pvalCols 

for(i in 1:length(pvalCols)){
method = substr(pvalCols[i],3,20)
print(method)
tmp =myQQ(method)
print(tmp)
}
