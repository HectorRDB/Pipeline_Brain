library(here)
library(stringr)
library(merger)
library(clusterExperiment)
library(mclust)
library(readr)

datasets <- c("SMARTer_cells_MOp", "SMARTer_nuclei_MOp",
              "10x_cells_MOp", "10x_nuclei_MOp")
# datasets <- c("SMARTer_cells_MOp",  "SMARTer_nuclei_MOp", )


for (dataset in datasets) {
  print(dataset)
  merger <- readRDS(here("data", type(dataset),
                         paste0(dataset, "_no_allen_mergers.rds")))
  
  print("...Initial consensus")
  clusters <- mergers_no_allen$initalMat
  r <- which(colnames(clusters) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    clusters[,"Rsec"] <- assignRsec(mergers_no_allen) 
  }
  clusters <- as.matrix(clusters) 
  cellsConsensus <- Consensus(clusMat = clusters,
                              large = (type != "Smart-Seq"))
  consensusInit <- cellsConsensus
  
  print("...Final consensus")
  clusters <- mergers_no_allen$currentMat
  r <- which(colnames(clusters) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    clusters[,"Rsec"] <- assignRsec(mergers_no_allen) 
  }
  clusters <- as.matrix(clusters) 
  
  cellsConsensus <- Consensus(clusMat = clusters,
                              large = (type != "Smart-Seq"))
  consensusFinal <- cellsConsensus$clustering
  
  print("...Intermediary consensus at 33.3%")
  midMat <- intermediateMat(merger = mergers_no_allen,
                            p = 1/3)
  r <- which(colnames(midMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    midMat[,"Rsec"] <- assignRsec(mergers_no_allen, p = 1/3)
  }
  midMat <- as.matrix(midMat)
  
  cellsConsensus <- Consensus(clusMat = midMat,
                              large = (type != "Smart-Seq"))
  consensusInt_33 <- cellsConsensus
  
  print("...Intermediary consensus at 66.7%")
  midMat <- intermediateMat(merger = mergers_no_allen,
                            p = 2/3)
  r <- which(colnames(midMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    midMat[,"Rsec"] <- assignRsec(mergers_no_allen, p = 2/3)
  }
  midMat <- as.matrix(midMat)
  
  cellsConsensus <- Consensus(clusMat = midMat,
                              large = (type != "Smart-Seq"))
  consensusInt_66 <- cellsConsensus
  
  print("...Intermediary consensus at 90%")
  midMat <- intermediateMat(merger = mergers_no_allen,
                            p = .9)
  r <- which(colnames(midMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    midMat[,"Rsec"] <- assignRsec(mergers_no_allen, p = .9)
  }
  midMat <- as.matrix(midMat)
  
  cellsConsensus <- Consensus(clusMat = midMat,
                              large = (type != "Smart-Seq"))
  consensusInt_90 <- cellsConsensus
  
  print("...Full matrix")
  names <- read_csv(here("data", type(dataset),
                         paste0(dataset, "_cluster.membership.csv")))
  clusters <- mergers_no_allen$initalMat
  r <- which(colnames(clusters) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    clusters[,"Rsec"] <- assignRsec(mergers_no_allen) 
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

