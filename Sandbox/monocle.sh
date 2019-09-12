#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/10x_nuclei_MOp_zinb.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/10x_nuclei_MOp_monocle_all.csv"
Rscript --verbose  monocle.R -l $loc -o $out > test.out 2>&1
