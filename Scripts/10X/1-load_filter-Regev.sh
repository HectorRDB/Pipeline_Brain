#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p RM
#SBATCH -t 12:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/Regev/count_matrix.csv"
out="/pylon5/ib5phhp/hectorrb/Regev/count_matrix_filt.rds"

module load intel/18.4
module load hdf5

Rscript --no-save --verbose  1-load_filter.R -l $loc -o $out -c 25 > 1b.out 2>&1
