#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
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

# Zinbwave
out3="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_norm.rds"
plot3="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp"
cluster="/scratch/users/singlecell/MiniAtlas/data/SMARTer_nuclei_MOp/cluster.annotation.csv"
MEMORYFILE="3b-memoryLogger.txt"

# while true; do free -h >> $MEMORYFILE; sleep 15; done & \
# Rscript --vanilla --verbose  3-zinb.R -n 10 -i 8 -f 20 -l $out2 -o $out3 -p $plot3 -c $cluster > 3b.out 2>&1
echo "Step 3"

# sc3
out4="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_sc3.rds"
#Rscript --vanilla --verbose  4-sc3.R -n 20 -l $out2 -o $out4 > 4b.out 2>&1
echo "Step 4"

# Seurat
out5="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_seurat.rds"
# Rscript --vanilla --verbose  5-seurat.R -l $out2 -o $out5 > 5b.out 2>&1
echo "Step 5"

# Monocle
out6="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_monocle.rds"
Rscript --vanilla --verbose 6-monocle.R -l $out3 -o $out6> 6b.out 2>&1
echo "Step 6"

# RSEC
MEMORYFILE="7b-memoryLogger.txt"
out7="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
# while true; do free -h >> $MEMORYFILE; sleep 15; done & \
# Rscript --vanilla --verbose  7-RSEC.R -n 30 -l $out3 -o $out7 > 7b.out 2>&1
echo "Step 7"

# ConsensusClustering
MEMORYFILE="8ab-memoryLogger.txt"
loc8="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --vanilla --verbose  7-ConsensusClustering.R -n 30 -l $loc8 -o $out -p $plot> 8ab.out 2>&1
echo "Step 8"
