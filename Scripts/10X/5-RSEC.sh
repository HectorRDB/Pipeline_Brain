#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
MEMORYFILE="5a_memoryLogger.txt"
loc="/pylon5/ib5phhp/hectorrb/out/10x_cells_MOp_zinbWs.rds"
out="/pylon5/ib5phhp/hectorrb/out/10x_cells_MOp_RSEC.rds"

module load hdf5

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  5-RSEC.R -n 10 -l $loc -o $out> 5a.out 2>&1
