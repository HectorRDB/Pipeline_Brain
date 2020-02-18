#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -t 7-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

# Create the MEMORYFILE file
timestamp=$(date +"%Y%m%d-%H%M%S")
basename=sc3_10x-v3-cells_${timestamp}
MEMORYFILE=${basename}.txt
# Add the first few lines for analytic purposes.
# The name variable can also be defined gloably by modifyinh the .bashrc file
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 10 LM 500GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

loc="/pylon5/ib5phnp/hectorrb/ProcessedData/10x_v3_cells_MOp_filt.rds"
out="/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_v3_cells_MOp_SC3.csv"

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & Rscript \
  --no-save --verbose  3-sc3.R -n 1 -l $loc -o $out > 3d.out 2>&1

logStorage=/pylon5/ib5phnp/hectorrb/logs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
