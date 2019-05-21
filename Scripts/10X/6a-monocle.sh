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

timestamp=$(date +"%m-%d-%H:%M")
MEMORYFILE=${timestamp}_6a-monocle_n-1.txt

echo $NAME > $MEMORYFILE
echo 1 "LM" 1000GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

# Example of a script you could run
# This can be more than one line
source activate monocle_env
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
  R CMD BATCH  6a-monocle.R 6a.out

cp $MEMORYFILE ${logStorage}/$MEMORYFILE
