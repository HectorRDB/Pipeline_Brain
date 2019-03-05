#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
MEMORYFILE="6a_memoryLogger.txt"
loc="/pylon5/ib5phhp/hectorrb/out/10x_cells_MOp"

module load hdf5

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  6-ConsensusClustering.R -n 10 -l $loc -o $loc > 6a.out 2>&1
