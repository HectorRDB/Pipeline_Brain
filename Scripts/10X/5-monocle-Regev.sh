#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

# This is not run on bridges but on the scf cluster
loc="/scratch/users/singlecell/MiniAtlas/data/rds/Regev_zinb.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/Regev_Monocle.csv"
Rscript --verbose  5-monocle.R -l $loc -o $out > 5-Regev.out 2>&1
