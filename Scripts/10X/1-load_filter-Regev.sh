#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p RM
#SBATCH -t 12:00:00
#SBATCH --nodes=1

module load intel/18.4
module load hdf5
module load gcc

timestamp=$(date +"%Y%m%d-%H%M%S")
basename=1-load-filter_Regev_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE


loc="/pylon5/ib5phhp/hectorrb/Regev/count_matrix.csv"
out="/pylon5/ib5phhp/hectorrb/Regev/count_matrix_filt.rds"

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
  Rscript --no-save --verbose  1-load_filter.R -l $loc -o $out -c 25 > 1-Regev.out 2>&1


logStorage=/pylon5/ib5phnp/hectorrb/logs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
