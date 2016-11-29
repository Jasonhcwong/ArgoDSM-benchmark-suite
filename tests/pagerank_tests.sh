#!/bin/bash
i=0
n=4
v=2000000
max=2
while [ $i -lt $max ]
do
	./../scripts/job.sh pagerank $n $v | sbatch
	echo "Job started for pagerank_${v} on $n nodes"
	true $((n=n*2))
	true $((i=i+1))
done
