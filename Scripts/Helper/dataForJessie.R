library(here)
library(stringr)
library(clusterExperiment)
library(mclust)
library(readr)

# datasets <- "SMARTer_cells_MOp SMARTer_nuclei_MOp 10x_cells_MOp 10x_nuclei_MOp"
source(here("Report", "helper.R"))
datasets <- c("SMARTer_cells_MOp",  "SMARTer_nuclei_MOp")


for (dataset in datasets) {
  print(dataset)
  merger <- readRDS(here("data", type(dataset),
                         paste0(dataset, "_no_allen_mergers.rds")))
  
  print("...Initial consensus")
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = as.matrix(merger$initalMat[, -2]),
                  clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 100)
  )
  consensusInit <- cellsConsensus$clustering
  
  print("...Final consensus")
  currentMat <- merger$currentMat
  currentMat[, "Rsec"] <- assignRsec(merger)
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = currentMat, 
                  clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 100)
  )
  consensusFinal <- cellsConsensus$clustering
  
  print("...Intermediary consensus at 33.3%")
  stopMatrix_33 <- intermediateMat(merger = merger, p = 1/3)
  stopMatrix_33[, "Rsec"] <- assignRsec(merger, p = 1/3)
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = merger$currentMat, clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 100)
  )
  consensusInt_33 <- cellsConsensus$clustering
  
  print("...Intermediary consensus at 66.7%")
  stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
  stopMatrix_66[, "Rsec"] <- assignRsec(merger, p = 2/3)
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = merger$currentMat, clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 100)
  )
  consensusInt_66 <- cellsConsensus$clustering
  
  print("...Intermediary consensus at 90%")
  stopMatrix_90 <- intermediateMat(merger = merger)
  stopMatrix_90[, "Rsec"] <- assignRsec(merger, p = .9)
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = merger$currentMat, clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 100)
  )
  consensusInt_90 <- cellsConsensus$clustering
  
  print("...Full matrix")
  names <- read_csv(here("data", type(dataset),
                         paste0(dataset, "_cluster.membership.csv")))
  currentMat <- merger$currentMat
  j <- which(colnames(merger$initalMat) == "RsecT")
  mat <- cbind(names$X1,
               merger$initalMat[, -j], consensusInit,
               stopMatrix_33, consensusInt_33,
               stopMatrix_66, consensusInt_66,
               stopMatrix_90, consensusInt_90,
               currentMat, consensusFinal)
  
  colnames(mat) <- c("cells",
    paste(c("sc3", "RSEC", "Monocle", "Seurat", "Consensus"), "Initial",
          sep = "-"),
    paste(c("sc3", "RSEC", "Monocle", "Seurat", "Consensus"), "33",
          sep = "-"),
    paste(c("sc3", "RSEC", "Monocle", "Seurat", "Consensus"), "66",
          sep = "-"),
    paste(c("sc3", "RSEC", "Monocle", "Seurat", "Consensus"), "90",
          sep = "-"),
    paste(c("sc3", "RSEC", "Monocle", "Seurat", "Consensus"), "Final",
          sep = "-")
    )
  
  write_csv(x = mat,
            path = here("data", "ForJessie", paste0(dataset, ".csv")))
}

