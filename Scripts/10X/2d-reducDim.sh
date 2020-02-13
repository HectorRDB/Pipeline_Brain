#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

module load gcc/8.2.0

# Create the MEMORYFILE file
timestamp=$(date +"%Y%m%d-%H%M%S")
basename=zinb_10x-v3-cells_${timestamp}
MEMORYFILE=${basename}.txt
# Add the first few lines for analytic purposes.
# The name variable can also be defined gloably by modifyinh the .bashrc file
NAME=Hector
echo $NAME > $MEMORYFILE
# Replace with your own variables. This is cpus-per-tasks partition mem
echo 20 LM 500GB >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

loc="/pylon5/ib5phnp/hectorrb/ProcessedData/10x_v3_cells_MOp_filt.rds"
out="/pylon5/ib5phnp/hectorrb/ProcessedData/10x_v3_cells_MOp_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/EDA/10x_v3_cells_Mop_tsne"
cluster="/pylon5/ib5phhp/hectorrb/10x_v3_cells_MOp/cluster.annotation.csv"

while true; do free -h >> $MEMORYFILE; sleep 30; done & Rscript --no-save --verbose\
  2-reducDim.R -l $loc -o $out -p $plot -n 8 -d 5 -c $cluster > 2d.out 2>&1

logStorage=/pylon5/ib5phnp/hectorrb/logs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
