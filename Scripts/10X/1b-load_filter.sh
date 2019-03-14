#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -t 24:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/10x_nuclei_MOp/"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_filt.rds"

module load intel/18.4
module load hdf5

Rscript --no-save --verbose  1-load_filter.R -l $loc -o $out -c 30 > 1a.out 2>&1
