#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J blastn
#SBATCH -n 5
#SBATCH -o /var/scratch/atangomo/analysis/teststdout.txt
#SBATCH -e /var/scratch/atangomo/analysis/teststderr.txt

# load the blast module
module load bcftools/1.8

bfctools query 
