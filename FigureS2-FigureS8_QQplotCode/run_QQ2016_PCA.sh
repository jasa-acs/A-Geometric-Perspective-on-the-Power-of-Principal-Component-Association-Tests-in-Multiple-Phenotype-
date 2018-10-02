modName=M8
mod1=8
sbatch -J M8 -t 6-12:00 --mem=20000  -p hsph   -o ./log/${modName}.log  R CMD BATCH --no-save "--args mod1=${mod1}" ./QQplot2016_PCA.R  ./log/${modName}.Rout




