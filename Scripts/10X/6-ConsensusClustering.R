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

# Load Monocle clustering results
print("Loading Monocle")
Monocle <- readRDS(paste0(loc, "_monocle2.rds"))

# Load sc3  and allen clustering results
print("Loading sc3")
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
print(colnames(colData(sc3)))
allen <- colData(sc3)[, "allenClusters"]
Monocle <- Monocle[rownames(colData(sc3))]
sc3 <- colData(sc3)[, "sc3_100_clusters"]


# Load Seurat
print("Loading Seurat")
seurat <- readRDS(paste0(loc, "_seurat.rds"))
seurat <- seurat[, "1.2,30"]

# Get the final clustering labels
  clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "seurat" = seurat)


# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
cat("Finished Consensus Merge\n")
saveRDS(object = mergers, file = paste0(output, "_mergers.rds"))