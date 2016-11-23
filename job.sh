#!/bin/bash -l
algorithm=$1
nodes=$2
argonodes=$(expr ${nodes} \* 2)
vertices=$3
threads=$(expr ${argonodes} \* 8)

echo "#!/bin/bash -l
#SBATCH -A g2016027
#SBATCH -p node -N ${nodes} -n ${threads} 
#SBATCH -t 5:00
#SBATCH -J ${algorithm}_${vertices} 
module load gcc/5.2.0
module load openmpi/1.8.8
cd /proj/g2016027/tim/project/${algorithm}
mpirun -v --map-by ppr:2:node --mca mpi_leave_pinned 1 --mca btl openib,self,sm -n ${argonodes} ./${algorithm} ${threads} testfiles/${algorithm}_${vertices}.txt output/${algorithm}_${vertices}.txt"


