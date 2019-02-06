#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out_norm="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_norm.rds"
out_rd="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_zinbWs.rds"
Rscript --vanilla --verbose  3-zinb.R -n 8 -l $loc -o $out_norm -r $out_rd> 3a.out 2>&1
