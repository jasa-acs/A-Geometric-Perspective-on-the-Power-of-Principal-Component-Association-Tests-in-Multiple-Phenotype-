mod1=8 
sbatch  -t 5-24:00 --mem=40000 -J M8 -p hsph  -o ./log/M8.log  R CMD BATCH --no-save "--args mod1=${mod1} " ./M.R ./log/M8.Rout


#mod1=10 
#sbatch  -t 5-24:00 --mem=40000 -J M10 -p hsph  -o ./log/M10.log  R CMD BATCH --no-save "--args mod1=${mod1} " ./M.R ./log/M10.Rout

