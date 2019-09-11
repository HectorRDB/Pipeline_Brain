#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_filt.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp_SC3.csv"
Rscript --verbose  4-sc3.R -n 20 -l $loc -o $out > 4b.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp_SC3.csv"
Rscript --verbose  4-sc3.R -n 20 -l $loc -o $out > 4a.out 2>&1
