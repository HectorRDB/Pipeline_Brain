#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_filt.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMerge/SMARTer_nuclei_MOp_singleSC3.csv"
Rscript --verbose  singleSC3.R -n 20 -l $loc -o $out> nuclei.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMerge/SMARTer_cells_MOp_singleSC3.csv"
Rscript --verbose  singleSC3.R -n 20 -l $loc -o $out> cells.out 2>&1
