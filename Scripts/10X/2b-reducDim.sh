#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p LM
#SBATCH --mem=1000GB
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0 
loc="/pylon5/ib5phhp/hectorrb/out/10x_cells_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/out/10x_cells_MOp_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/Exploration/10x_cells_Mop_tsne"
MEMORYFILE="2_zinb_memoryLogger.txt"

while true; do free -h >> $MEMORYFILE; sleep 30; done & Rscript --no-save --verbose\
  2b-reducDim.R -l $loc -o $out -p $plot -n 20 -d 1 > 2a.out 2>&1
