#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
Rscript --verbose  2-filtering.R -l $loc -o $out > 2a.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_filt.rds"
Rscript --verbose  2-filtering.R -l $loc -o $out > 2b.out 2>&1
