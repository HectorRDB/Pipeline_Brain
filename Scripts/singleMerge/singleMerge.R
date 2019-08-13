library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(tidyverse)
library(clusterExperiment)

loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp"
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
Rsec <- assignUnassigned(Rsec, clusterLabel = "allAssigned")
Rsec <- makeDendrogram(Rsec, whichCluster = "allAssigned")
cutoffs <- seq(from = .05, to = 1, by = .05)
names(cutoffs) <- cutoffs
res_nuclei <- map_df(cutoffs,
                     function(i){
                 Rsec2 <- mergeClusters(Rsec,
                                        mergeMethod = "adjP",
                                        plotInfo = "adjP",
                                        cutoff = i,
                                        clusterLabel = "Clusters",
                                        plot = F,
                                        DEMethod = "limma")  
                 
                 return(Rsec2@clusterMatrix[,"Clusters"])  
  })

loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_cells_MOp"
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
Rsec <- assignUnassigned(Rsec, clusterLabel = "allAssigned")
Rsec <- makeDendrogram(Rsec, whichCluster = "allAssigned")
res_cells <- map_df(cutoffs,
               function(i){
                 Rsec2 <- mergeClusters(Rsec,
                                        mergeMethod = "adjP",
                                        plotInfo = "adjP",
                                        cutoff = i,
                                        clusterLabel = "Clusters",
                                        plot = F,
                                        DEMethod = "limma")  
                 
              return(Rsec2@clusterMatrix[,"Clusters"])
               })

write_csv(res_cells, 
          path = "../../data/Smart-Seq/SMARTer_cells_MOp_Rsec_single_merge.csv")
write_csv(res_nuclei,
          path = "../../data/Smart-Seq/SMARTer_nuclei_MOp_Rsec_single_merge.csv")
