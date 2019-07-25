#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p LM
#SBATCH --mem=1000GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc
timestamp=$(date +"%m-%d-%H:%M")
basename=nuclei-cells-ARI-merging_allen-no-allen_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 10 "LM" 1000GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

echo "Nuclei dataset"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp"
out="/home/hectorrb/Pipeline_Brain/data/10X/10x_nuclei_MOp"
# while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
#       Rscript --no-save --verbose  7-ConsensusClustering.R -n 10 -l $loc \
#       -o $out -a FALSE > 7bb.out 2>&1

echo "Cells dataset"
loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp"
out="/home/hectorrb/Pipeline_Brain/data/10X/10x_cells_MOp"
while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
      Rscript --no-save --verbose  7-ConsensusClustering.R -n 10 -l $loc \
      -o $out -a FALSE > 7ab.out 2>&1

logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
