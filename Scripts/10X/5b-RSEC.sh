#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH -p LM
#SBATCH --mem=1500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

MEMORYFILE="5b-memoryLogger.txt"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/10X/10x_nuclei_Mop_clusterMany.pdf"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_RSEC.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  6-RSEC.R -n 15 -l $loc -o $out -p $plot > 5b.out 2>&1
