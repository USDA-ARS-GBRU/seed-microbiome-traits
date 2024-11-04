#!/bin/bash
#SBATCH --job-name=mvmiss
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --partition=medium
#SBATCH --time=7-00:00:00

cd /home/quentin.read/GitHub/seed-microbiome-traits
module load r
Rscript2 brm_mv_model_testing_effectrange.R

