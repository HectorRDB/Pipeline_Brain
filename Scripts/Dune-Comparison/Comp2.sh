#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_cells_MOp"
rsec="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/SMARTer_cells_MOp_large2"
Rscript --verbose  Dune.R -n 20 -l $loc -o $out -S "1.7.30" -C "15" -m "k_15" -r $rsec > cells2.out 2>&1

loc="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleMethod/SMARTer_nuclei_MOp"
rsec="/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_Rsec.rds"
out="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/SMARTer_nuclei_MOp_large2"
Rscript --verbose  Dune.R -n 20 -l $loc -o $out -S "1.7.30" -C "15" -m "k_15" -r $rsec > nuclei2.out 2>&1
