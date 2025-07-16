#!/bin/sh

#SBATCH --account=psych
#SBATCH --job-name=supervispr
#SBATCH --output=logs/supervisor-%j
#SBATCH --time=11:59:00

# loop through subjects in elfk_mem_sublist.txt and run fmriprep on them
participants=`cat elfk_mem_sublist.txt`

for ppt in ${participants[@]}
do 
	echo $ppt
	sbatch ./scripts_ginsburg/run_fmriprep.sh $ppt
	sleep 10m
done
