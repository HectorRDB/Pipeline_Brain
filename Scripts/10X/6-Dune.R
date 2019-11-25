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
colnames(sc3) <- str_remove(colnames(sc3), "^X")
Names <- sc3$cells
sc3 <- sc3[, sc3_p] %>% as.numeric()

# Load Seurat clustering results
print("Loading Seurat")
seurat <- read.csv(paste0(loc, "_Seurat.csv"))[, -1]
colnames(seurat) <- str_remove(colnames(seurat), "^X")
ggsave(filename = paste0(opt$p, "_Seurat_ARI.png"),
       plot = plotARIs(seurat %>% select(-cells)))
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
print("Loading Monocle")
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
ggsave(filename = paste0(opt$p, "_Monocle_ARI.png"),
       plot = plotARIs(Monocle %>% select(-cells)))
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Get the final clustering labels
Names <- as.character(Names)
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = seurat)
rownames(clusMat) <- Names

# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
merger <- Dune(clusMat = clusMat, nCores = opt$n)
cat("Finished Consensus Merge\n")
saveRDS(object = merger, file = paste0(output, "_mergers.rds"))

# Save the matrix with all the consensus steps ----
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

write_csv(x = as.data.frame(mat), path = paste0(output, ".csv"))