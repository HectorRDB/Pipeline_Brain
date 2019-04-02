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
                  proportion = 2/3, minSize = 20)
  )
  consensusInit <- cellsConsensus$clustering
  
  print("...Final consensus")
  currentMat <- merger$currentMat
  currentMat[, 2] <- assignRsec(merger)
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = currentMat, 
                  clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 20)
  )
  consensusFinal <- cellsConsensus$clustering
  
  print("...Intermediary consensus")
  stopMatrix <- intermediateMat(merger = merger)
  stopMatrix[, "Rsec"] <- assignRsec(merger, p = .9)
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = merger$currentMat, clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 20)
  )
  consensusInt <- cellsConsensus$clustering
  
  print("...Full matrix")
  names <- read_csv(here("data", type(dataset),
                         paste0(dataset, "_cluster.membership.csv")))
  mat <- cbind(names$X1,
               merger$initalMat[,-2], consensusInit,
               stopMatrix, consensusInt,
               currentMat, consensusFinal)
  
  colnames(mat) <- c("cells",
    paste(c("sc3", "RSEC", "Seurat", "Consensus"), "Initial", sep = "-"),
    paste(c("sc3", "RSEC", "Seurat", "Consensus"), "Intermediate", sep = "-"),
    paste(c("sc3", "RSEC", "Seurat", "Consensus"), "Final", sep = "-")
    )
  
  write_csv(x = mat,
            path = here("data", "ForJessie", paste0(dataset, ".csv")))
}

