#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

R CMD BATCH --no-save 4b-sc3.R 4b-sc3.out
