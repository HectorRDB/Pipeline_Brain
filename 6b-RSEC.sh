#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

Rscript --vanilla --verbose  6-RSEC.R -d "SMARTer_nuclei_MOp/" > 6b.out 2>&1
