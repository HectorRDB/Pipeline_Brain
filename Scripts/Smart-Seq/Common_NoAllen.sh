#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --nodes=1

# Nuclei
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_nuclei_MOp_no_allen"

Rscript --vanilla --verbose  8-ConsensusClustering.R -n 30 -a FALSE -l $loc -o $out > 8bb.out 2>&1
echo "Nuclei dataset"

# Cell
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Smart-Seq/SMARTer_cells_MOp_no_allen"

Rscript --vanilla --verbose  8-ConsensusClustering.R -n 30 -a FALSE -l $loc -o $out > 8ba.out 2>&1
echo "Cell dataset"
