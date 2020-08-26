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
sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
colnames(sc3) <- str_remove(colnames(sc3), "^X") %>% str_replace("\\.", "-")
ggsave(filename = paste0(opt$p, "_SC3_NMI.png"),
       plot = plotNMIs(sc3 %>% select(-cells), values = FALSE,
                       numericalLabels = TRUE))
Names <- sc3$cells
sc3 <- sc3[,"0"]

# Load Seurat clustering results
seurat <- read.csv(paste0(loc, "_Seurat.csv"))[, -1]
colnames(seurat) <- str_remove(colnames(seurat), "^X")
ggsave(filename = paste0(opt$p, "_Seurat_NMI.png"),
       plot = plotNMIs(seurat %>% select(-cells), values = FALSE))
seurat_p <- "1.2.50"
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
ggsave(filename = paste0(opt$p, "_Monocle_NMI.png"),
       plot = plotNMIs(Monocle %>% select(-cells), values = FALSE))
monocle_p <- "k_45"
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Get the final clustering labels
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = seurat)
rownames(clusMat) <- Names  

# Do the consensus clustering ----
BPPARAM <- BiocParallel::MulticoreParam(opt$n)
print(paste0("Number of cores: ", opt$n))
print(system.time(
  merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, unclustered = -1,
                 verbose = TRUE, parallel = TRUE, metric = "NMI")
))
cat("Finished Consensus Merge\n")
saveRDS(object = merger, file = paste0(output, "_mergersNMI.rds"))

# Save the matrix with all the consensus steps ----
Names <- as.character(Names)
chars <- c("sc3", "Monocle", "Seurat")

levels <- seq(from = 0, to = 1, by = .05)
stopMatrix <- lapply(levels, function(p){
  print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
  mat <- intermediateMat(merger = merger, p = p) %>%
    as.matrix()
  mat <- mat[Names, ]
  return(mat)
}) %>%
  do.call('cbind', args = .)

colnames(stopMatrix) <- lapply(levels, function(p){
  i <- as.character(round(100 * p))
  if (nchar(i) == 1) {
    i <- paste0("0", i)
  }
  return(paste(chars, i, sep = "-"))
}) %>% unlist()
print("...Full matrix")
mat <- cbind(as.character(Names), stopMatrix)

colnames(mat)[1] <- "cells"


write_csv(x = as.data.frame(mat), path = paste0(output, "_NMI.csv"))