#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc

timestamp=$(date +"%m-%d-%H:%M")
MEMORYFILE=${timestamp}_dataForJessie_4-datasets.txt

echo $NAME > $MEMORYFILE
echo 1 "LM" 500GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

# Example of a script you could run
# This can be more than one line
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
  R CMD BATCH  dataForJessie.R dataForJessie.out

cp $MEMORYFILE ${logStorage}/$MEMORYFILE
