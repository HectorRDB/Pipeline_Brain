#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

R CMD BATCH --no-save 4a-sc3C.R 4a-sc3C.out
R CMD BATCH --no-save 4a-sc3N.R 4a-sc3N.out
