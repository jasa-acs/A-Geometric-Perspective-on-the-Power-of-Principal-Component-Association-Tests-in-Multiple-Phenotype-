#modName=M9_1PCO 
mod1=9 
for chunkid  in {1..100}
do
sbatch  -t 1-24:00 --mem=3000 -J M9_${chunkid} -p general  -o ./log/M9_${chunkid}.log  R CMD BATCH --no-save "--args mod1=${mod1}  chunkid=${chunkid}" ./pcaqpco.R ./log/M9_${chunkid}.Rout
done 
#modName=M10_1PCO 
mod1=10 
for chunkid in {1..100}
do 
sbatch  -t 1-24:00 --mem=3000 -J M10_${chunkid} -p general  -o ./log/M10_${chunkid}.log  R CMD BATCH --no-save "--args mod1=${mod1}  chunkid=${chunkid}" ./pcaqpco.R ./log/M10_${chunkid}.Rout
done
