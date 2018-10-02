
################################################################
# Read_in all the data 
################################################################
## SBP 
if(SBP_in == TRUE){
dat.SBP = read.table(SBP,header=T)
dat.SBP$Z  = dat.SBP$beta/dat.SBP$se
dat.SBP = dat.SBP[is.na(dat.SBP$Z)==FALSE,] ## keep non-mis_ing rows 
dat.SBP = dat.SBP[is.nan(dat.SBP$Z)==FALSE,]
p = dim(dat.SBP)[2]
post = ".SBP"
names(dat.SBP)[2:p] = paste0(names(dat.SBP)[2:p],rep(post,p-1))
}

## DBP 
if(DBP_in == TRUE){
dat.DBP  = read.table(DBP,header=T)
dat.DBP$Z  = dat.DBP$beta/dat.DBP$se
dat.DBP = dat.DBP[is.na(dat.DBP$Z)==FALSE,] ## keep non-mis_ing rows 
dat.DBP = dat.DBP[is.nan(dat.DBP$Z)==FALSE,] ## remove not an number rows 
p = dim(dat.DBP)[2]
post = ".DBP"
names(dat.DBP)[2:p] = paste0(names(dat.DBP)[2:p],rep(post,p-1))
}


## BMI 
if(BMI_in == TRUE){
dat.BMI = read.table(BMI,header=T)
dat.BMI$Z  = dat.BMI$beta/dat.BMI$se
dat.BMI = dat.BMI[is.na(dat.BMI$Z)==FALSE,] ## keep non-mis_ing rows 
dat.BMI = dat.BMI[is.nan(dat.BMI$Z)==FALSE,]
p = dim(dat.BMI)[2]
post = ".BMI"
names(dat.BMI)[2:p] = paste0(names(dat.BMI)[2:p],rep(post,p-1))
}
## WHR 
if(WHR_in == TRUE){
dat.WHR = read.table(WHR,header=T)
dat.WHR$Z  = dat.WHR$beta/dat.WHR$se
dat.WHR = dat.WHR[is.na(dat.WHR$Z)==FALSE,] ## keep non-mis_ing rows 
dat.WHR = dat.WHR[is.nan(dat.WHR$Z)==FALSE,]
p = dim(dat.WHR)[2]
post = ".WHR"
names(dat.WHR)[2:p] = paste0(names(dat.WHR)[2:p],rep(post,p-1))
}

## fGLU 
if(fGLU_in == TRUE){
dat.fGLU = read.table(fGLU,header=T)
dat.fGLU$Z  = dat.fGLU$beta/dat.fGLU$se
dat.fGLU = dat.fGLU[is.na(dat.fGLU$Z)==FALSE,] ## keep non-mis_ing rows 
dat.fGLU = dat.fGLU[is.nan(dat.fGLU$Z)==FALSE,]
p = dim(dat.fGLU)[2]
post = ".fGLU"
names(dat.fGLU)[2:p] = paste0(names(dat.fGLU)[2:p],rep(post,p-1))
}

## fINSULIN 
if(fINSULIN_in == TRUE){
dat.fINSULIN = read.table(fINSULIN,header=T)
dat.fINSULIN$Z  = dat.fINSULIN$beta/dat.fINSULIN$se
dat.fINSULIN = dat.fINSULIN[is.na(dat.fINSULIN$Z)==FALSE,] ## keep non-mis_ing rows 
dat.fINSULIN = dat.fINSULIN[is.nan(dat.fINSULIN$Z)==FALSE,]
p = dim(dat.fINSULIN)[2]
post = ".fINSULIN"
names(dat.fINSULIN)[2:p] = paste0(names(dat.fINSULIN)[2:p],rep(post,p-1))
}

## LDL 
if(LDL_in == TRUE){
dat.LDL = read.table(LDL,header=T)
dat.LDL$Z  = dat.LDL$beta/dat.LDL$se
dat.LDL = dat.LDL[is.na(dat.LDL$Z)==FALSE,] ## keep non-mis_ing rows 
dat.LDL = dat.LDL[is.nan(dat.LDL$Z)==FALSE,]
p = dim(dat.LDL)[2]
post = ".LDL"
names(dat.LDL)[2:p] = paste0(names(dat.LDL)[2:p],rep(post,p-1))
}

## TG 
if(TG_in == TRUE){
dat.TG = read.table(TG,header=T)
dat.TG$Z  = dat.TG$beta/dat.TG$se
dat.TG = dat.TG[is.na(dat.TG$Z)==FALSE,] ## keep non-mis_ing rows 
dat.TG = dat.TG[is.nan(dat.TG$Z)==FALSE,]
p = dim(dat.TG)[2]
post = ".TG"
names(dat.TG)[2:p] = paste0(names(dat.TG)[2:p],rep(post,p-1))
}

### HDL 
if(HDL_in == TRUE){
dat.HDL = read.table(HDL,header=T)
dat.HDL$Z  = dat.HDL$beta/dat.HDL$se
dat.HDL = dat.HDL[is.na(dat.HDL$Z)==FALSE,] ## keep non-mis_ing rows 
dat.HDL = dat.HDL[is.nan(dat.HDL$Z)==FALSE,]
p = dim(dat.HDL)[2]
post = ".HDL"
names(dat.HDL)[2:p] = paste0(names(dat.HDL)[2:p],rep(post,p-1))
}

## TC
if(TC_in == TRUE){
dat.TC = read.table(TC,header=T)
dat.TC$Z  = dat.TC$beta/dat.TC$se
dat.TC = dat.TC[is.na(dat.TC$Z)==FALSE,] ## keep non-mis_ing rows 
dat.TC = dat.TC[is.nan(dat.TC$Z)==FALSE,]
p = dim(dat.TC)[2]
post = ".TC"
names(dat.TC)[2:p] = paste0(names(dat.TC)[2:p],rep(post,p-1))
}


