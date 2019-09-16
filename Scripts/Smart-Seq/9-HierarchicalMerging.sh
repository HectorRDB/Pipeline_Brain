#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

# Cells
loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp"
rsec="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/SMARTer_cells_MOp_hierarchical.csv"
Rscript --verbose 9-HierarchicalMerging.R -l $loc -o $out -r $rsec > 9a.out 2>&1

# nuclei
loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp"
rsec="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/SMARTer_nuclei_MOp_hierarchical.csv"
Rscript --verbose 9-HierarchicalMerging.R -l $loc -o $out -r $rsec > 9b.out 2>&1
