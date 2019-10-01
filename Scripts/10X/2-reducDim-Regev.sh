#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -p LM
#SBATCH --mem=800GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

timestamp=$(date +"%Y%m%d-%H%M%S")
basename=2-zinbWave_Regev_${timestamp}
MEMORYFILE=${basename}.txt

echo $NAME > $MEMORYFILE
echo 1 "RM" >> $MEMORYFILE
TIMELAPSES=30
echo $TIMELAPSES >> $MEMORYFILE

module load gcc/8.2.0
loc="/pylon5/ib5phhp/hectorrb/Regev/count_matrix_filt.rds"
out="/pylon5/ib5phhp/hectorrb/Regev/count_matrix_norm.rds"
plot="/home/hectorrb/Pipeline_Brain/Figures/EDA/Regev_tsne"
cluster="/pylon5/ib5phhp/hectorrb/10x_nuclei_MOp/cluster.annotation.csv"

while true; do free -h >> $MEMORYFILE; sleep $TIMELAPSES; done & \
  2-reducDim.R -l $loc -o $out -p $plot -n 20 -d 5 -c $cluster > Regev.out 2>&1


logStorage=/pylon5/ib5phhp/shared/improved-happiness/xsede/xsedelogs
cp $MEMORYFILE ${logStorage}/$MEMORYFILE
