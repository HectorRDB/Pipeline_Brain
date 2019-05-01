#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
MEMORYFILE="6b-memoryLogger.txt"
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_zinbWs.rds"
out="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_monocle.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose 6-monocle.R -l $loc -o $out> 6b.out 2>&1
