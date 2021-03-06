#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p LM
#SBATCH --mem=1500GB
#SBATCH -t 24:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_filt.rds"
out_norm="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_norm.rds"
out_rd="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_zinbWs"
MEMORYFILE="2a_memoryLogger.txt"
MEMORYSUMMARY="2a_memorySummary.txt"

module load cuda/8.0 pytorch/0.1.5 python3/intel_3.6.3

while true; do free -h >> $MEMORYFILE; sleep 30; done & \
./2-scVI.py -n 5 -l $loc -o $out_norm -r $out_rd> 2a.out 2>&1
