#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
MEMORYFILE="7a-memoryLogger.txt"
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_zinbWs.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  7-RSEC.R -n 10 -l $loc -o $out> 7a.out 2>&1
