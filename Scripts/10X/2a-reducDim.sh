#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p LM
#SBATCH --mem=1500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/Exploration/10x_cells_Mop_tsne"
cluster="/pylon5/ib5phhp/hectorrb/10x_cells_MOp/cluster.annotation.csv"
MEMORYFILE="2_zinb_memoryLogger.txt"

while true; do free -h >> $MEMORYFILE; sleep 30; done & Rscript --no-save --verbose\
  2-reducDim.R -l $loc -o $out -p $plot -n 8 -d 5 -c $cluster > 2a.out 2>&1
