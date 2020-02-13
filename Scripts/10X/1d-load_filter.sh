#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM
#SBATCH -t 24:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/10x_v3_cells_MOp"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_v3_cells_MOp_filt.rds"

timestamp=$(date +"%Y%m%d-%H%M%S")
basename=1-load-filter_10x-v3-cells_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

module load intel/18.4
module load hdf5
module load gcc

Rscript --no-save --verbose  1-load_filter.R -l $loc -o $out -c 30 > 1d.out 2>&1
