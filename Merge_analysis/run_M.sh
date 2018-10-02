modName=M8_1 
mod1=8
mod2=1 
SBP_in=TRUE 
BMI_in=TRUE 
WHR_in=TRUE 
fINSULIN_in=TRUE 
fGLU_in=TRUE 
LDL_in=TRUE
TG_in=TRUE 
HDL_in=TRUE 
sbatch  -t 0-12:00 --mem=10000 -p serial_requeue  -o ./log/${modName}.log  R CMD BATCH --no-save "--args mod1=${mod1} mod2=${mod2}  SBP_in=${SBP_in}  BMI_in=${BMI_in} WHR_in=${WHR_in} fINSULIN_in=${fINSULIN_in} fGLU_in=${fGLU_in} LDL_in=${LDL_in} TG_in=${TG_in} HDL_in=${HDL_in}" ./M.R ./log/${modName}.Rout

