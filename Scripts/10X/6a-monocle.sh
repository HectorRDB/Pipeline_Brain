#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
MEMORYFILE="6ab-memoryLogger.txt"
loc="/scratch/users/singlecell/MiniAtlas/data/rds/10x_cells_MOp"
out="/scratch/users/singlecell/MiniAtlas/data/rds/10x_cells_MOp_monocle2.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose 6-monocle.R -l $loc -o $out> 6ab.out 2>&1
