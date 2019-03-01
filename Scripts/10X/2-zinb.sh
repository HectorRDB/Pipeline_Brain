#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_filt.rds"
out_rd="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_zinbWs"
MEMORYFILE="2a_memoryLogger.txt"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  2-zinb.R -n 4 -l $loc -r $out_rd -t TRUE > 2a.out 2>&1
