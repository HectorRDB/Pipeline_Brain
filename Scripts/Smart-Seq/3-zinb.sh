#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_zinbWs.rds"
MEMORYFILE="3a_memoryLogger.txt"
MEMORYSUMMARY="3a_memorySummary.txt"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  3-zinb.R -n 4 -l $loc -o $out> 3a.out 2>&1
