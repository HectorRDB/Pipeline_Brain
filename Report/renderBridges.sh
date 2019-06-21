#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1

# Create the MEMORYFILE file
timestamp=$(date +"%Y%m%d-%H%M%S")
# Do not use underscore in your script name or your parameters
basename=render-files_SMART-Seq_${timestamp}
MEMORYFILE=${basename}.txt

# Add the first few lines for analytic purposes.
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 1 RM 128GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

module load gcc
module load pandoc
datasets="SMARTer_cells_MOp SMARTer_nuclei_MOp"

for dataset in $datasets
do
  echo $dataset
  while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
    Rscript --vanilla --verbose render.R -d $dataset -l FALSE \
      > ${basename}.out 2>&1
done

# Note that if your script fail, you will need to execute this line manually.
logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
