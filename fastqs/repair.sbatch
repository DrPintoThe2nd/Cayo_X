#!/bin/bash
#SBATCH -N 1  				# number of nodes
#SBATCH -n 8  				# number of "tasks" (default: allocates 1 core per task)
#SBATCH -t 0-8:00:00   		# time in d-hh:mm:ss
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

#repair.sh in=GTEX-11GSP-0004-SM-6WBTQ.fq.gz out1=GTEX-11GSP-0004-SM-6WBTQ_R1.fq.gz out2=GTEX-11GSP-0004-SM-6WBTQ_R2.fq.gz outs=GTEX-11GSP-0004-SM-6WBTQ_singletons.fq.gz repair;

#repair.sh in=GTEX-131YS-0003-SM-7DRPE.fq.gz out1=GTEX-131YS-0003-SM-7DRPE_R1.fq.gz out2=GTEX-131YS-0003-SM-7DRPE_R2.fq.gz outs=GTEX-131YS-0003-SM-7DRPE_singletons.fq.gz repair;

#repair.sh in=GTEX-14DAQ-0004-SM-7DROH.fq.gz out1=GTEX-14DAQ-0004-SM-7DROH_R1.fq.gz out2=GTEX-14DAQ-0004-SM-7DROH_R2.fq.gz outs=GTEX-14DAQ-0004-SM-7DROH_singletons.fq.gz repair;

repair.sh in=GTEX-15ER7-0002-SM-7DLFM.fq.gz out1=GTEX-15ER7-0002-SM-7DLFM_R1.fq.gz out2=GTEX-15ER7-0002-SM-7DLFM_R2.fq.gz outs=GTEX-15ER7-0002-SM-7DLFM_singletons.fq.gz repair;

