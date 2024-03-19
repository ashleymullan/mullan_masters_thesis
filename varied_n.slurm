#!/bin/bash
#SBATCH --job-name="varied_n"          # Name that appears in queue
#SBATCH --partition=small               # Resource group (small/medium/large)
#SBATCH --account=lotspeichGrp          # Research group
#SBATCH --nodes=1                       # Number of Nodes
#SBATCH --ntasks-per-node=1             # Number of CPUS
#SBATCH --mem=16G                       # Requested memory
#SBATCH --time=00-06:00:00              # Job duration in DD-HH:MM:SS
#SBATCH --output="varied_n-%j.out"     # Slurm stdout, %j is the job number

hostname

module load R/4.2.1

# Change to directory where .R script exists
cd /deac/sta/lotspeichGrp/mullae22/thesis

# Execute add.m script
R CMD BATCH varied_n_full_sims.R

slurm_mem_report -g
