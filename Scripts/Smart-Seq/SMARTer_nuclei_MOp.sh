#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

# Load Data
loc="/scratch/users/singlecell/MiniAtlas/data/SMARTer_nuclei_MOp/"
out1="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp.rds"
# Rscript --vanilla --verbose  1-loadData.R -l $loc -o $out1 > 1b.out 2>&1
echo "Step 1"

# Filter data
out2="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_filt.rds"
# Rscript --vanilla --verbose  2-filtering.R -l $out1 -o $out2 -c 120 -r 60 > 2b.out 2>&1
echo "Step 2"

# sc3
out4="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_sc3.rds"
#Rscript --vanilla --verbose  4-sc3.R -n 20 -l $out2 -o $out4 > 4b.out 2>&1
echo "Step 4"

# Seurat
out5="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_seurat.rds"
# Rscript --vanilla --verbose  5-seurat.R -l $out2 -o $out5 > 5b.out 2>&1
echo "Step 5"

# Zinbwave
out3="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_norm.rds"
MEMORYFILE="3b_memoryLogger.txt"
MEMORYSUMMARY="3b_memorySummary.txt"

# while true; do free -h >> $MEMORYFILE; sleep 15; done & \
# Rscript --vanilla --verbose  3-zinb.R -n 10 -l $out2 -o $out3> 3b.out 2>&1
echo "Step 3"

# RSEC
MEMORYFILE="6b_memoryLogger.txt"
out6="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"

# while true; do free -h >> $MEMORYFILE; sleep 15; done & \
# Rscript --vanilla --verbose  6-RSEC.R -n 20 -l $out -o $out6 -k 3 > 6b.out 2>&1
echo "Step 6"

# ConsensusClustering
MEMORYFILE="7b_memoryLogger.txt"
loc7="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  7-ConsensusClustering.R -n 20 -l $loc7 -o $loc7 -p $plot> 7b.out 2>&1
echo "Step 7"
