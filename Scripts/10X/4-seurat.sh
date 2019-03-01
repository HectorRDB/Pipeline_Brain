#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_seurat.rds"
Rscript --vanilla --verbose  4-seurat.R -l $loc -o $out> 4a.out 2>&1
