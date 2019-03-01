#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/10x_cells_MOp/"
out="/pylon5/ib5phhp/hectorrb/rds/10x_cells_MOp_filt.rds"
Rscript --vanilla --verbose  1-load_filter.R -l $loc -o $out> 1a.out 2>&1
