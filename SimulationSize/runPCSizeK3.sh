
modName=unstr1v2
K=3
rho=0.3
alpha=0.05
simNum=1e7
unstr=TRUE
ncore=50
sbatch  -t 6-12:00 --mem-per-cpu=4000 -p serial_requeue -n 50  -J ${modName} -o ./log/${modName}.log  R CMD BATCH --no-save "--args K=${K} ncore=${ncore} rho=${rho} alpha=${alpha} simNum=${simNum} unstr=${unstr}" PCSizeK3.R  ./log/${modName}.Rout


modName=unstr2v2
K=3
rho=0.3
alpha=0.01
simNum=1e7
unstr=TRUE
ncore=50
sbatch  -t 6-12:00 --mem-per-cpu=4000 -p serial_requeue -n 50  -J ${modName} -o ./log/${modName}.log  R CMD BATCH --no-save "--args K=${K} ncore=${ncore} rho=${rho} alpha=${alpha} simNum=${simNum} unstr=${unstr}" PCSizeK3.R  ./log/${modName}.Rout



modName=unstr3v2
K=3
rho=0.3
alpha=0.001
simNum=1e7
unstr=TRUE
ncore=50
sbatch  -t 6-12:00 --mem-per-cpu=4000 -p serial_requeue -n 50  -J ${modName} -o ./log/${modName}.log  R CMD BATCH --no-save "--args K=${K} ncore=${ncore} rho=${rho} alpha=${alpha} simNum=${simNum} unstr=${unstr}" PCSizeK3.R  ./log/${modName}.Rout

modName=unstr4v2
K=3
rho=0.3
alpha=0.0001
simNum=1e7
unstr=TRUE
ncore=50
sbatch  -t 6-12:00 --mem-per-cpu=4000 -p serial_requeue -n 50  -J ${modName} -o ./log/${modName}.log  R CMD BATCH --no-save "--args K=${K} ncore=${ncore} rho=${rho} alpha=${alpha} simNum=${simNum} unstr=${unstr}" PCSizeK3.R  ./log/${modName}.Rout

modName=unstr5v2
K=3
rho=0.3
alpha=0.00001
simNum=1e7
unstr=TRUE
ncore=50
sbatch  -t 6-12:00 --mem-per-cpu=4000 -p serial_requeue -n 50  -J ${modName} -o ./log/${modName}.log  R CMD BATCH --no-save "--args K=${K} ncore=${ncore} rho=${rho} alpha=${alpha} simNum=${simNum} unstr=${unstr}" PCSizeK3.R  ./log/${modName}.Rout

