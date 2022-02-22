#!/bin/bash
#SBATCH -N 1  				# number of nodes
#SBATCH -n 8  				# number of "tasks" (default: allocates 1 core per task)
#SBATCH -t 0-12:00:00   		# time in d-hh:mm:ss
#SBATCH -p serial	       		# partition 
#SBATCH -q normal       		# QOS
#SBATCH -o slurm.%j.out 		# file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err 		# file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL 		# Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user=bpinto2@asu.edu  	# Mail-to address
#SBATCH --job-name=umap 		#Job name
#SBATCH --export=NONE

module purge

unset PERL5LIB
export PERL5LIB=$PERL5LIB:/home/bpinto2/miniconda3/envs/snakemake/lib/site_perl

source activate Cayo

bash unmap_cram.sh

