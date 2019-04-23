#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

MEMORYFILE="5b-memoryLogger.txt"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/10X/10x_nuclei_Mop_clusterMany.pdf"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_RSEC.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --no-save --verbose  5-RSEC.R -n 5 -l $loc -o $out -p $plot > 5b.out 2>&1
