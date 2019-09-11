suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1, type = "integer",
              help = "Number of cores to use"
  ),
  make_option(c("-p", "--plot"),
              action = "store", default = 1, type = "integer",
              help = "Where to store the plots"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at ", loc, "\n")
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
  cat("The output will be stored at ", output, "\n")
}

library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(tidyverse)
library(Dune)

# Load Data ----
# Load sc3 clustering results
sc3 <- read.csv(paste0(loc, "_SC3.csv"))
ggsave(filename = paste0(opt$p, "_monocle_ARI.png"),
       plot = clusterMatToAri(sc3 %>% select(-cells)))
Names <- sc3$cells
sc3 <- sc3[,"sc3_0_clusters"]

# Load Seurat clustering results
seurat <- read.csv(paste0(loc, "_seurat.csv"))
ggsave(filename = paste0(opt$p, "_seurat_ARI.png"),
       plot = clusterMatToAri(seurat %>% select(-cells)))
seurat_p <- "1.2,50"
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))
ggsave(filename = paste0(opt$p, "_monocle_ARI.png"),
       plot = clusterMatToAri(Monocle %>% select(-cells)))
monocle_p <- "k_45"
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Load RSEC clustering results
Rsec <- read.csv(paste0(loc, "Rsec.csv"))
RsecT <- Rsec$RsecT
Rsec <- Rsec$Rsec

# Get the final clustering labels
clusMat <- data.frame("sc3" = sc3, "Rsec" = Rsec, "Monocle" = Monocle,
                      "seurat" = seurat)
rownames(clusMat) <- Names  

# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
print(system.time(
  merger <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
))
cat("Finished Consensus Merge\n")
merger$initalMat <- cbind(mergers$initalMat, RsecT)
colnames(merger$initalMat)[5] <- "RsecT"
saveRDS(object = merger, file = paste0(output, "_mergers.rds"))

# Save the matrix with all the consensus steps ----
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

write_csv(x = as.data.frame(mat), path = paste0(output, ".csv"))