#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

Rscript --vanilla --verbose  5-Seurat.R -d "SMARTer_cells_MOp/" > 5a.out 2>&1
