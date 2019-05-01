#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
MEMORYFILE="7a_memoryLogger.txt"
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_cells_MOp_mergers.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  7-ConsensusClustering.R -n 10 -l $loc -o $loc > 7a.out 2>&1
