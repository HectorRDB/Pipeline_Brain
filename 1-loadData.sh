#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/SMARTer_cells_MOp/"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp.rds"
Rscript --vanilla --verbose  1-loadData.R -l $loc -o $out> 1a.out 2>&1
