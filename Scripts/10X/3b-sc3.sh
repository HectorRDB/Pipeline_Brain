#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p LM
#SBATCH --mem=800GB
#SBATCH -t 7-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

# Create the MEMORYFILE file
timestamp=$(date +"%Y%m%d-%H%M%S")
basename=sc3_10x-nuclei_${timestamp}
MEMORYFILE=${basename}.txt
# Add the first few lines for analytic purposes.
# The name variable can also be defined gloably by modifying the .bashrc file
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 1 LM 800GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_sc3.rds"

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & Rscript \
  --no-save --verbose  3-sc3.R -n 1 -l $loc -o $out> ${basename}.out 2>&1

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
