#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

MEMORYFILE=replicability.txt

while true; do free -h >> $MEMORYFILE; sleep 30 ; done & \
  R CMD BATCH metaneighbor_analysis.R metaneighbor_analysis.Rout

while true; do free -h >> $MEMORYFILE; sleep 30 ; done & \
  R CMD BATCH visualization.R visualization.Rout