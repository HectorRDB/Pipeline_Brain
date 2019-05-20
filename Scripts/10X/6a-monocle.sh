#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p LM
#SBATCH --mem=1000GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_norm.rds"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_monocle.rds"
module unload intel
module load anaconda3

# Create the MEMORYFILE file
timestamp=$(date +"%m-%d-%H:%M")
MEMORYFILE=${timestamp}_6a-monocle_n-1.txt

# Add the first few lines for analytic purposes.
# The name variable can also be defined gloably by modifyinh the .bashrc file
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 1 "LM" 1000GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

# Example of a script you could run
# This can be more than one line
source activate monocle_env
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
  R CMD BATCH --vanilla --verbose 6-monocle.R 

# This can also be stored globally
logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
