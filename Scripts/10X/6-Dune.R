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
  make_option(c("-p", "--plot"),
              action = "store", default = NA, type = "character",
              help = "Where to save the plots"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1, type = "integer",
              help = "Number of cores to use"
  ),
  make_option(c("-S", "--SeuratParam"),
              action = "store", default = NA, type = "character",
              help = "Parameter to use for Seurat"
  ),
  make_option(c("-C", "--C3"),
              action = "store", default = NA, type = "character",
              help = "SC3 parameter"
  ),
  make_option(c("-m", "--monocle"),
              action = "store", default = NA, type = "character",
              help = "Monocle parameter"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc)
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
}

if (!is.na(opt$C)) {
  sc3_p <- opt$C
} else {
  stop("Missing C argument")
}

if (!is.na(opt$C)) {
  seurat_p <- opt$S
} else {
  stop("Missing S argument")
}

if (!is.na(opt$m)) {
  monocle_p <- opt$m
} else {
  stop("Missing m argument")
}

# Load Data----
library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(Dune)

# Load sc3 clustering results
print("Loading SC3")
sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
Names <- sc3$cells
sc3 <- sc3[, paste0("sc3_", sc3_p, "_clusters")] %>% as.numeric()

# Load Seurat clustering results
print("Loading Seurat")
seurat <- read.csv(paste0(loc, "_Seurat.csv"))[, -1]
colnames(seurat) <- str_remove(colnames(seurat), "^X")
ggsave(filename = paste0(opt$p, "_Seurat_ARI.png"),
       plot = clusterMatToAri(seurat %>% select(-cells)))
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
print("Loading Monocle")
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
ggsave(filename = paste0(opt$p, "_Monocle_ARI.png"),
       plot = clusterMatToAri(Monocle %>% select(-cells)))
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Get the final clustering labels
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = seurat)

# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
merger <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
cat("Finished Consensus Merge\n")
saveRDS(object = merger, file = paste0(output, "_mergers.rds"))

# Save the matrix with all the consensus steps ----
print("...Initial consensus")
initialMat <- as.matrix(merger$initalMat) 
cellsConsensus <- Consensus(clusMat = initialMat, large = TRUE)
consensusInit <- cellsConsensus

print("...Final consensus")
currentMat <- merger$currentMat
currentMat <- as.matrix(currentMat) 

cellsConsensus <- Consensus(clusMat = currentMat, large = TRUE)
consensusFinal <- cellsConsensus

print("...Intermediary consensus at 33.3%")
stopMatrix_33 <- intermediateMat(merger = merger,
                                 p = 1/3)
stopMatrix_33 <- as.matrix(stopMatrix_33)

cellsConsensus <- Consensus(clusMat = stopMatrix_33, large = TRUE)
consensusInt_33 <- cellsConsensus

print("...Intermediary consensus at 66.7%")
stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
stopMatrix_66 <- as.matrix(stopMatrix_66)

cellsConsensus <- Consensus(clusMat = stopMatrix_66, large = TRUE)
consensusInt_66 <- cellsConsensus

print("...Intermediary consensus at 90%")
stopMatrix_90 <- intermediateMat(merger = merger, p = .9)
stopMatrix_90 <- as.matrix(stopMatrix_90)

cellsConsensus <- Consensus(clusMat = stopMatrix_90, large = TRUE)
consensusInt_90 <- cellsConsensus

print("...Full matrix")
mat <- cbind(as.character(Names),
             initialMat, consensusInit,
             stopMatrix_33, consensusInt_33,
             stopMatrix_66, consensusInt_66,
             stopMatrix_90, consensusInt_90,
             currentMat, consensusFinal)
chars <- c("sc3", "Monocle", "Seurat", "Consensus")

colnames(mat) <- c("cells",
                   paste(chars, "Initial", sep = "-"), paste(chars, "33", sep = "-"),
                   paste(chars, "66", sep = "-"), paste(chars, "90", sep = "-"),
                   paste(chars, "Final", sep = "-")
)

write_csv(x = as.data.frame(mat), path = paste0(output, ".csv"))