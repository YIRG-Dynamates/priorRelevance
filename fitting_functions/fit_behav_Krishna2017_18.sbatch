#!/usr/bin/env bash

# == SLURM directives ==
#SBATCH --job-name=behav_Krishna2017_18

# docs https://docs.vbc.ac.at/books/scientific-computing/page/cbe-slurm-users-advanced-guide

#SBATCH --nodes=2
#SBATCH --ntasks=30
#SBATCH --ntasks-per-node=15
#SBATCH --cpus-per-task=2
#SBATCH --time=0-07:59:00         
#SBATCH --qos=short
#SBATCH --mem=80G
 
#SBATCH --output=behav_Krishna2017_18_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.meijer@services.vbc.ac.at

# rm -rf ~/.matlab
rm ~/java.log.*

# Load Matlab
ml build-env/f2021
ml matlab/2021b_update_2

# Launch job, `srun` knows how to do this
srun -l  matlab -nodesktop -nosplash -r "fit_behav_Krishna2017_18"

# Schedule next job if needed
sbatch fit_behav_Krishna2017_19.sbatch
