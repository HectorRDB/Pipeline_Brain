#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

R CMD BATCH --no-save 5a-seurat_cluster.R 5a-seurat_cluster.out
