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
# loc <- "data/singleMethod/SMARTer_cells_MOp"
# opt <- list(p = "Figures/Smart-Seq/SMARTer_cells_MOp", n = 1)
# output <- "data/Dune/SMARTer_cells_MOp"

# Load Data ----
# Load sc3 clustering results
sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
colnames(sc3) <- str_remove(colnames(sc3), "^X") %>% str_replace("\\.", "-")
ggsave(filename = paste0(opt$p, "_SC3_ARI.png"),
       plot = plotARIs(sc3 %>% select(-cells), values = FALSE,
                       numericalLabels = TRUE))
Names <- sc3$cells
sc3 <- sc3[,"0"]

# Load Seurat clustering results
seurat <- read.csv(paste0(loc, "_Seurat.csv"))[, -1]
colnames(seurat) <- str_remove(colnames(seurat), "^X")
ggsave(filename = paste0(opt$p, "_Seurat_ARI.png"),
       plot = plotARIs(seurat %>% select(-cells), values = FALSE))
seurat_p <- "1.2.50"
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
ggsave(filename = paste0(opt$p, "_Monocle_ARI.png"),
       plot = plotARIs(Monocle %>% select(-cells), values = FALSE))
monocle_p <- "k_45"
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Get the final clustering labels
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = seurat)
rownames(clusMat) <- Names  

# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
print(system.time(
  merger <- Dune(clusMat = clusMat, nCores = opt$n, unclustered = -1,
                 verbose = TRUE)
))
cat("Finished Consensus Merge\n")
saveRDS(object = merger, file = paste0(output, "_mergers.rds"))

# Save the matrix with all the consensus steps ----
print("...Initial consensus")
initialMat <- as.matrix(merger$initalMat) 
cellsConsensus <- Consensus(clusMat = initialMat, large = FALSE)
consensusInit <- cellsConsensus

print("...Final consensus")
currentMat <- merger$currentMat
currentMat <- as.matrix(currentMat) 

cellsConsensus <- Consensus(clusMat = currentMat, large = FALSE)
consensusFinal <- cellsConsensus

print("...Intermediary consensus at 33.3%")
stopMatrix_33 <- intermediateMat(merger = merger,
                                 p = 1/3)
stopMatrix_33 <- as.matrix(stopMatrix_33)

cellsConsensus <- Consensus(clusMat = stopMatrix_33, large = FALSE)
consensusInt_33 <- cellsConsensus

print("...Intermediary consensus at 66.7%")
stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
stopMatrix_66 <- as.matrix(stopMatrix_66)

cellsConsensus <- Consensus(clusMat = stopMatrix_66, large = FALSE)
consensusInt_66 <- cellsConsensus

print("...Intermediary consensus at 90%")
stopMatrix_90 <- intermediateMat(merger = merger, p = .9)
stopMatrix_90 <- as.matrix(stopMatrix_90)

cellsConsensus <- Consensus(clusMat = stopMatrix_90, large = FALSE)
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