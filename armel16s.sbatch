#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J Mothur_Job
#SBATCH -n 5
#SBATCH -w taurus
#SBATCH -o /var/scratch/atangomo/analysis/stdout.txt
#SBATCH -e /var/scratch/atangomo/analysis/stderr.txt


#####SBATCH -M tangomoarmel25@gmail.com

#module load openmpi/1.8.2
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load mothur/1.39.5

cd /var/scratch/atangomo/analysis/
mothur 1.run_mothur.batch
