#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1

MEMORYFILE="8a-memoryLogger.txt"
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_cells_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_cells_MOp"
while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --verbose  8-ConsensusClustering.R -n 10 -l $loc -o $out -p $plot > 7a.out 2>&1

MEMORYFILE="8b-memoryLogger.txt"
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_nuclei_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp"
while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --verbose  8-ConsensusClustering.R -n 10 -l $loc -o $out -p $plot > 7b.out 2>&1
