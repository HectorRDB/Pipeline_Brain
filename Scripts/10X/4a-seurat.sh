#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p LM
#SBATCH --mem=1000GB
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1

loc="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_filt.rds"
out="/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_seurat.rds"
module load intel/18.4
module load hdf5


Rscript --no-save --verbose  4-seurat.R -l $loc -o $out> 4a.out 2>&1
