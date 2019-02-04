#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_sc3.rds"
Rscript --vanilla --verbose  4-sc3.R -n 8 -l $loc -o $out> 4a.out 2>&1
