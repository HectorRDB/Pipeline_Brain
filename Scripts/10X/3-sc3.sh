#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_sc3.rds"
Rscript --vanilla --verbose  3-sc3.R -n 8 -l $loc -o $out> 3a.out 2>&1
