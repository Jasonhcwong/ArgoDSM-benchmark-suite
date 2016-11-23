#!/bin/bash -l
algorithm=$1
nodes=$2
argonodes=$(expr ${nodes} \* 2)
vertices=$3
threads=$(expr ${argonodes} \* 8)

echo "#!/bin/bash -l
#SBATCH -p cluster -N ${nodes} -n ${threads} 
#SBATCH -t 5:00
#SBATCH -J ${algorithm}_${vertices} 
cd ../${algorithm}
mpirun -v --map-by ppr:2:node --mca mpi_leave_pinned 1 --mca btl openib,self,sm -n ${argonodes} ./${algorithm} ${threads} testfiles/${algorithm}_${vertices}.txt output/${algorithm}_${vertices}.txt"


