#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p RM
#SBATCH -t 24:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

# Create the MEMORYFILE file
timestamp=$(date +"%Y%m%d-%H%M%S")
basename=rep_all_${timestamp}
MEMORYFILE=${basename}.txt
# Add the first few lines for analytic purposes.
# The name variable can also be defined gloably by modifyinh the .bashrc file
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 20 LM 500GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

while true; do free -h >> $MEMORYFILE; sleep 30 ; done & \
  R CMD BATCH --no-save 02-metaneighbor_analysis.R metaneighbor_analysis.Rout

while true; do free -h >> $MEMORYFILE; sleep 30 ; done & \
  R CMD BATCH --no-save 03-visualization.R visualization.Rout

logStorage=/pylon5/ib5phnp/hectorrb/logs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
