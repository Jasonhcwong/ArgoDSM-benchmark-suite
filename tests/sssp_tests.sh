#!/bin/bash
i=0
n=1
v=1024
max=4
while [ $i -lt $max ]
do
	./../scripts/sssp_job.sh sssp $n $v | sbatch
	echo "Job started for sssp_${v} on $n nodes"
	true $((n=n*2))
	true $((i=i+1))
done
