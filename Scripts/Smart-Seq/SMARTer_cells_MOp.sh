#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

# Load Data
loc="/scratch/users/singlecell/MiniAtlas/data/SMARTer_cells_MOp/"
out1="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp.rds"
# Rscript --vanilla --verbose  1-loadData.R -l $loc -o $out1 > 1a.out 2>&1
echo "Step 1"

# Filter data
out2="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_filt.rds"
# Rscript --vanilla --verbose  2-filtering.R -l $out1 -o $out2> 2a.out 2>&1
echo "Step 2"

# sc3
out4="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_sc3.rds"
#Rscript --vanilla --verbose  4-sc3.R -n 20 -l $out2 -o $out4 > 4a.out 2>&1
echo "Step 4"

# Seurat
out5="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_seurat.rds"
# Rscript --vanilla --verbose  5-seurat.R -l $out2 -o $out5 > 5a.out 2>&1
echo "Step 5"

# Zinbwave
out3="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_norm.rds"
plot3="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_cells_MOp"
cluster="/scratch/users/singlecell/MiniAtlas/data/SMARTer_cells_MOp/cluster.annotation.csv"
MEMORYFILE="3a_memoryLogger.txt"

# while true; do free -h >> $MEMORYFILE; sleep 15; done & \
# Rscript --vanilla --verbose  3-zinb.R -n 6 -l $out2 -o $out3 -p $plot3 -c $cluster > 3a.out 2>&1
echo "Step 3"

# RSEC
MEMORYFILE="6a_memoryLogger.txt"
out6="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"

# while true; do free -h >> $MEMORYFILE; sleep 15; done & \
# Rscript --vanilla --verbose  6-RSEC.R -n 20 -l $out -o $out6 -k 3 > 6a.out 2>&1
echo "Step 6"

# ConsensusClustering
MEMORYFILE="7a_memoryLogger.txt"
loc7="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_cells_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_cells_MOp"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  7-ConsensusClustering.R -n 20 -l $loc7 -o $out -p $plot> 7a.out 2>&1
echo "Step 7"
