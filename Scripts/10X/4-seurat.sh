#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load intel/18.4
module load hdf5
module load gcc
timestamp=$(date +"%m-%d-%H:%M")
basename=seurat_cells-nuclei_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

echo "cells dataset"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_filt.rds"
out="/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_cells_MOp_Seurat.csv"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --no-save --verbose  4-seurat.R -l $loc -o $out> 4a.out 2>&1

echo "Nuclei dataset"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_filt.rds"
out="/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_nuclei_MOp_Seurat.csv"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --no-save --verbose  4-seurat.R -l $loc -o $out> 4b.out 2>&1

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
