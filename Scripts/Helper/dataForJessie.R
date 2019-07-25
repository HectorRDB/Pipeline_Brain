library(here)
library(stringr)
library(merger)
library(clusterExperiment)
library(mclust)
library(readr)

datasets <- c("SMARTer_cells_MOp", "SMARTer_nuclei_MOp",
              "10x_cells_MOp", "10x_nuclei_MOp")
# datasets <- c("SMARTer_cells_MOp",  "SMARTer_nuclei_MOp", )
types <- function(dataset) {
  if (str_detect(dataset, "SMART")) return("Smart-Seq")
  if (str_detect(dataset, "10x")) return("10X")
  stop("Type unknown")
}

for (dataset in datasets) {
  print(dataset)
  type <- types(dataset)
  merger <- readRDS(here("data", types(dataset),
                         paste0(dataset, "_no_allen_mergers.rds")))
  
  print("...Initial consensus")
  clusters <- merger$initalMat
  r <- which(colnames(clusters) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    clusters[,"Rsec"] <- assignRsec(merger) 
  }
  clusters <- as.matrix(clusters) 
  cellsConsensus <- Consensus(clusMat = clusters,
                              large = (type != "Smart-Seq"))
  consensusInit <- cellsConsensus
  
  print("...Final consensus")
  clusters <- merger$currentMat
  r <- which(colnames(clusters) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    clusters[,"Rsec"] <- assignRsec(merger) 
  }
  clusters <- as.matrix(clusters) 
  
  cellsConsensus <- Consensus(clusMat = clusters,
                              large = (type != "Smart-Seq"))
  consensusFinal <- cellsConsensus$clustering
  
  print("...Intermediary consensus at 33.3%")
  midMat <- intermediateMat(merger = merger,
                            p = 1/3)
  r <- which(colnames(midMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    midMat[,"Rsec"] <- assignRsec(merger, p = 1/3)
  }
  midMat <- as.matrix(midMat)
  
  cellsConsensus <- Consensus(clusMat = midMat,
                              large = (type != "Smart-Seq"))
  consensusInt_33 <- cellsConsensus
  
  print("...Intermediary consensus at 66.7%")
  midMat <- intermediateMat(merger = merger,
                            p = 2/3)
  r <- which(colnames(midMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    midMat[,"Rsec"] <- assignRsec(merger, p = 2/3)
  }
  midMat <- as.matrix(midMat)
  
  cellsConsensus <- Consensus(clusMat = midMat,
                              large = (type != "Smart-Seq"))
  consensusInt_66 <- cellsConsensus
  
  print("...Intermediary consensus at 90%")
  midMat <- intermediateMat(merger = merger,
                            p = .9)
  r <- which(colnames(midMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    midMat[,"Rsec"] <- assignRsec(merger, p = .9)
  }
  midMat <- as.matrix(midMat)
  
  cellsConsensus <- Consensus(clusMat = midMat,
                              large = (type != "Smart-Seq"))
  consensusInt_90 <- cellsConsensus
  
  print("...Full matrix")
  names <- read_csv(here("data", types(dataset),
                         paste0(dataset, "_cluster.membership.csv")))
  clusters <- merger$initalMat
  r <- which(colnames(clusters) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    clusters[,"Rsec"] <- assignRsec(merger) 
  }
  clusters <- as.matrix(clusters) 
  mat <- cbind(names$X1,
               clusters, consensusInit,
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

