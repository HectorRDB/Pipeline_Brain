#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/SMARTer_cells_MOp"
Rscript --verbose  8-ConsensusClustering.R -n 20 -l $loc -o $out -S "1.7,30" -C "10" -m "k_20" > cells.out 2>&1

loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/SMARTer_nuclei_MOp"
Rscript --verbose  8-ConsensusClustering.R -n 20 -l $loc -o $out -S "1.7,30" -C "10" -m "k_20" > nuclei.out 2>&1
