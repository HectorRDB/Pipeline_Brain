#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM
#SBATCH -t 24:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phnp/hectorrb/10x_v3_nuclei_MOp/"
out="/pylon5/ib5phnp/hectorrb/ProcessedData/10x_v3_nuclei_MOp_filt2.rds"

timestamp=$(date +"%Y%m%d-%H%M%S")
basename=1-load-filter_10x-v3-nuclei_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

module load intel/18.4
module load hdf5
module load gcc

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
    Rscript --no-save --verbose  1-load_filter.R -l $loc -o $out \
    -c 10 > 1c.out 2>&1


logStorage=/pylon5/ib5phnp/hectorrb/logs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
