#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

# Cell
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_no_allen"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_cells_MOp_no_allen"

Rscript --vanilla --verbose  7-ConsensusClustering.R -n 32 -a FALSE\
        -l $loc -o $loc -p $plot > 7ab.out 2>&1
echo "Cell dataset"

# Nuclei
loc="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_no_allen"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp_no_allen"

Rscript --vanilla --verbose  7-ConsensusClustering.R -n 32 -a FALSE\
        -l $loc -o $loc -p $plot> 7bb.out 2>&1
echo "Nuclei dataset"
