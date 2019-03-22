#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p LM
#SBATCH --mem=1500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_sc3.rds"
MEMORYFILE="3a-memoryLogger.txt"

while true; do free -h >> $MEMORYFILE; sleep 30; done & Rscript \
  --no-save --verbose  3-sc3.R -n 1 -l $loc -o $out> 3a.out 2>&1
