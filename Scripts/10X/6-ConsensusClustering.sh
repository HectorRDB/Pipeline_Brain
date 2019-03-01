#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
MEMORYFILE="6a_memoryLogger.txt"
loc="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp"
plots="~/Pipeline_Brain/Figures/10x/10x_cells_MOp"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  6-ConsensusClustering.R -n 10 -l $loc -o $loc -p $plots > 6a.out 2>&1
