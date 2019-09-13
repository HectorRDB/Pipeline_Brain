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
library(dplyr)
library(tidyr)
library(stringr)
library(Dune)

# Load sc3
print("Loading sc3")
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
Names <- colnames(sc3)
sc3 <- colData(sc3)[, sc3_p] %>% as.numeric()

# Load Seurat clustering results
seurat <- read.csv(paste0(loc, "_seurat.csv"))
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Get the final clustering labels
  clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "seurat" = seurat)


# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
cat("Finished Consensus Merge\n")
saveRDS(object = mergers, file = paste0(output, "_mergers.rds"))