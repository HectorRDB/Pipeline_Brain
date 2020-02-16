#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

# This is not run on bridges but on the scf cluster
# For the 10x v3 we can use the smart-seq script
Pipeline="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain"

loc="/scratch/users/singlecell/MiniAtlas/data/rds/10x_v3_nuclei_MOp_norm.rds"
out=${Pipeline}/singleMethod/10x_v3_nuclei_MOp_Monocle.csv
Rscript --verbose ${Pipeline}/Smart-Seq/6-monocle.R -l $loc -o $out > 5c.out 2>&1

loc="/scratch/users/singlecell/MiniAtlas/data/rds/10x_v3_cells_MOp_zinb.rds"
out=${Pipeline}/singleMethod/10x_v3_cells_MOp_Monocle.csv
Rscript --verbose ${Pipeline}/Smart-Seq/6-monocle.R -l $loc -o $out > 5d.out 2>&1

