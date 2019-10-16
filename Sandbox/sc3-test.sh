#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1

module load intel/18.4
module load hdf5
module load gcc
timestamp=$(date +"%Y%m%d-%H%M%S")
basename=sc3-explo_10x-cells_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      R CMD BATCH sc3-test.R sc3-test.Rout

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
