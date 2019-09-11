#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_filt.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_zinbWs.rds"
Rscript --verbose  3-zinb.R -n 20 -l $loc -o $out > 3b.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_zinbWs.rds"
Rscript --verbose  3-zinb.R -n 20 -l $loc -o $out > 3a.out 2>&1
