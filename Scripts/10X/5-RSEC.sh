#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p LM
#SBATCH --mem=1500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

MEMORYFILE="5a-memoryLogger.txt"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/Exploration/10x_cells_Mop_tsne"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_RSEC.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  6-RSEC.R -n 10 -l $loc -o $out -p $plot > 5a.out 2>&1
