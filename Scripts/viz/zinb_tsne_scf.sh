#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_zinbWs.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_cells_MOp_tnse.csv"
Rscript --vanilla --verbose zinb_tsne.R -l $loc -o $out -K 20 > cells.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_zinbWs.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_nuclei_MOp_tnse.csv"
Rscript --vanilla --verbose zinb_tsne.R -l $loc -o $out -K 11 > nuclei.out 2>&1
