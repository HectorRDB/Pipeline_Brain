library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(tidyverse)
library(clusterExperiment)

loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))

for (i in seq(from = .1, to = 1, by = .05)) {
  Rsec2 <- mergeClusters(Rsec,
                        mergeMethod = "adjP",
                        plotInfo = "adjP",
                        cutoff = i,
                        clusterLabel = "Clusters",
                        plot = F,
                        DEMethod = "limma")  
  print(i)
  print(n_distinct(Rsec2@clusterMatrix[,"Clusters"]))
}