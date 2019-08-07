library(here)
library(stringr)
library(merger)
library(clusterExperiment)
library(mclust)
library(readr)

datasets <- c("SMARTer_cells_MOp", "SMARTer_nuclei_MOp",
              "10x_cells_MOp", "10x_nuclei_MOp")
datasets <- c("SMARTer_cells_MOp",  "SMARTer_nuclei_MOp")
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
  initialMat <- merger$initalMat
  r <- which(colnames(initialMat) == "RsecT")
  if (all.equal(integer(0) ,r) != TRUE) {
    initialMat <- initialMat[, -r]
  }
  initialMat <- as.matrix(initialMat) 
  cellsConsensus <- Consensus(clusMat = initialMat,
                              large = (type != "Smart-Seq"))
  consensusInit <- cellsConsensus
  
  print("...Final consensus")
  currentMat <- merger$currentMat
  currentMat <- as.matrix(currentMat) 
  
  cellsConsensus <- Consensus(clusMat = currentMat,
                              large = (type != "Smart-Seq"))
  consensusFinal <- cellsConsensus
  
  print("...Intermediary consensus at 33.3%")
  stopMatrix_33 <- intermediateMat(merger = merger,
                            p = 1/3)
  stopMatrix_33 <- as.matrix(stopMatrix_33)
  
  cellsConsensus <- Consensus(clusMat = stopMatrix_33,
                              large = (type != "Smart-Seq"))
  consensusInt_33 <- cellsConsensus
  
  print("...Intermediary consensus at 66.7%")
  stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
  stopMatrix_66 <- as.matrix(stopMatrix_66)
  
  cellsConsensus <- Consensus(clusMat = stopMatrix_66,
                              large = (type != "Smart-Seq"))
  consensusInt_66 <- cellsConsensus
  
  print("...Intermediary consensus at 90%")
  stopMatrix_90 <- intermediateMat(merger = merger, p = .9)
  stopMatrix_90 <- as.matrix(stopMatrix_90)
  
  cellsConsensus <- Consensus(clusMat = stopMatrix_90,
                              large = (type != "Smart-Seq"))
  consensusInt_90 <- cellsConsensus
  
  print("...Full matrix")
  names <- read_csv(here("data", types(dataset),
                         paste0(dataset, "_cluster.membership.csv")))
  mat <- cbind(names$X1,
               initialMat, consensusInit,
               stopMatrix_33, consensusInt_33,
               stopMatrix_66, consensusInt_66,
               stopMatrix_90, consensusInt_90,
               currentMat, consensusFinal)
  if (type == "Smart-Seq") {
    chars <- c("sc3", "RSEC", "Monocle", "Seurat", "Consensus")
  } else {
    chars <- c("sc3", "Monocle", "Seurat", "Consensus")
  }
  
  colnames(mat) <- c("cells",
    paste(chars, "Initial", sep = "-"), paste(chars, "33", sep = "-"),
    paste(chars, "66", sep = "-"), paste(chars, "90", sep = "-"),
    paste(chars, "Final", sep = "-")
    )
  
  write_csv(x = as.data.frame(mat),
            path = here("data", "ForJessie", paste0(dataset, ".csv")))
}

