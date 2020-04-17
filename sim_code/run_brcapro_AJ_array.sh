#!/bin/bash 

#SBATCH -J brcaproAJ # A single job name for the array 
#SBATCH -n 16 # cores
#SBATCH -N 1 # nodes
#SBATCH -p serial_requeue # Partition 
#SBATCH -t 0-0:40 # Running time 
#SBATCH --mem 6500 # Memory request 
#SBATCH -o brcapro.out # Standard output 
#SBATCH -e brcapro.err # Standard error
#SBATCH --mail-type=END 
#SBATCH --mail-user=zguan@g.harvard.edu

module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER

seed=${SLURM_ARRAY_TASK_ID}
argString="--args "$seed

R CMD BATCH --quiet --no-restore --no-save "$argString" run_brcapro_AJ.R
