#!/bin/bash
#SBATCH -N 1  				# number of nodes
#SBATCH -n 1  				# number of "tasks" (default: allocates 1 core per task)
#SBATCH -t 0-01:00:00   		# time in d-hh:mm:ss
#SBATCH -p htc	        		# partition 
#SBATCH -q normal       		# QOS
#SBATCH -o logs/slurm.%j.out 		# file to save job's STDOUT (%j = JobId)
#SBATCH -e logs/slurm.%j.err 		# file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL 		# Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user=bpinto2@asu.edu  	# Mail-to address
#SBATCH --job-name=SMCT 		#Job name
#SBATCH --export=NONE

module purge

unset PERL5LIB
export PERL5LIB=$PERL5LIB:/home/bpinto2/miniconda3/envs/snakemake/lib/site_perl

source activate Cayo

snakemake --use-conda -np -s Snakefile_SNPs.py --unlock;

snakemake --use-conda -j 400 --rerun-incomplete --latency-wait 60 -s Snakefile_asereadcounter.py --cluster "sbatch -q normal -p htc -t 0-00:45:00 -N 1 -n 1 -o logs/slurm.snakemake.%j.out"
