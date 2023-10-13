#!/bin/bash
#SBATCH --job-name=nnetseed
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --partition=short
#SBATCH --time=2-00:00:00

cd /home/quentin.read/GitHub/seed-microbiome-traits
module load r/4.3.0 python_3
Rscript2 ${scriptname}

