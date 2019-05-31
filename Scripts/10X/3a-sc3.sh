#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p LM
#SBATCH --mem=1500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

# Create the MEMORYFILE file
timestamp=$(date +"%Y%m%d-%H%M%S")
basename=sc3_10x-cells_${timestamp}
MEMORYFILE=${basename}.txt
# Add the first few lines for analytic purposes.
# The name variable can also be defined gloably by modifyinh the .bashrc file
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 1 LM 1500GB >> $MEMORYFILE
TIMELAPSES=15
echo $TIMELAPSES >> $MEMORYFILE

loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_sc3.rds"

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & Rscript \
  --no-save --verbose  3-sc3.R -n 1 -l $loc -o $out > ${basename}.out 2>&1

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
