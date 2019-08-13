#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --nodes=1

# Nuclei
R CMD BATCH singleMerge.R
