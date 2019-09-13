#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p RM
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1

module load intel/18.4
module load hdf5
module load gcc
timestamp=$(date +"%m-%d-%H:%M")
basename=10x-consensus_cells-nuclei_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 20 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

echo "cells dataset"
loc="/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_cells_MOp"
out="/home/hectorrb/Pipeline_Brain/data/Dune/10x_cells_MOp"
plot="/home/hectorrb/Pipeline_Brain/Figures/10X/10x_cells_MOp"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --verbose  6-Dune.R -n 20 -l $loc -o $out -S "1.2,50" -C "100" \
          -m "k_45" -p $plot > 6a.out 2>&1

echo "Nuclei dataset"
loc="/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_cells_MOp"
out="/home/hectorrb/Pipeline_Brain/data/Dune/10x_cells_MOp"
plot="/home/hectorrb/Pipeline_Brain/Figures/10X/10x_cells_MOp"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --verbose  6-Dune.R -n 20 -l $loc -o $out -S "1.2,50" -C "100" \
          -m "k_45" -p $plot > 6b.out 2>&1

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE