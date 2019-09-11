#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_zinbWs.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp_Rsec.csv"
sce="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_Rsec.rds"
Rscript --verbose  7-RSEC.R -n 10 -l $loc -o $out -s $sce > 7b.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_zinbWs.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp_Rsec.csv"
sce="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_Rsec.rds"
Rscript --verbose  7-RSEC.R -n 10 -l $loc -o $out -s $sce > 7a.out 2>&1
