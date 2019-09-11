#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_filt.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp_Seurat.csv"
Rscript --verbose  5-seurat.R -l $loc -o $out > 5b.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp_Seurat.csv"
Rscript --verbose  5-seurat.R -l $loc -o $out > 5a.out 2>&1
