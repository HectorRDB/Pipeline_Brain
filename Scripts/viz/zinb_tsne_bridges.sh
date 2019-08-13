#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM
#SBATCH -t 0-20:00:00
#SBATCH --nodes=1

module load gcc
timestamp=$(date +"%m-%d-%H:%M")
basename=zinb-tsne_cells-nuclei_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

echo "Nuclei dataset"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_zinb.rds"
out="/home/hectorrb/Pipeline_Brain/data/10X/10x_nuclei_MOp_tnse.csv"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --no-save --verbose  zinb_tsne.R -l $loc -o $out -K 20 > nuclei.out 2>&1

echo "Cells dataset"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_zinb.rds"
out="/home/hectorrb/Pipeline_Brain/data/10X/10x_cells_MOp_tnse.csv"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --no-save --verbose  zinb_tsne.R -l $loc -o $out -K 30 > cells.out 2>&1

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
