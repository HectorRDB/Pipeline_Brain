#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Dune/SMARTer_cells_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_cells_MOp"
Rscript --verbose  8-Dune.R -n 20 -l $loc -o $out -p $plot > 8a.out 2>&1

loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Dune/SMARTer_nuclei_MOp"
plot="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Smart-Seq/SMARTer_nuclei_MOp"
Rscript --verbose  8-Dune.R -n 20 -l $loc -o $out -p $plot > 8b.out 2>&1
