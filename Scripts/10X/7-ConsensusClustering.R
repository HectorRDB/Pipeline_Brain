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
  make_option(c("-a", "--allen"),
              action = "store", default = T, type = "logical",
              help = "Wether to use allen or not"
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

library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(merger)

# Load Data and clean seurat ----
# Load sc3  and allen clustering results
print("Loading sc3")
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
allen <- colData(sc3)[, "allenClusters"]
sc3 <- colData(sc3)[, "sc3_100_clusters"]
print(sum(is.na(sc3)))
print("Allen")
print(sum(is.na(allen)))

# Load Monocle clustering results
print("Loading Monocle")
Monocle <- readRDS(paste0(loc, "_monocle2.rds"))
Monocle <- Monocle[names(sc3)]
print(sum(is.na(Monocle)))

# Load Seurat
print("Loading Seurat")
seurat <- readRDS(paste0(loc, "_seurat.rds"))
seurat <- seurat[, "1.2,30"]
print(sum(is.na(seurat)))

# Load Allen
# Get the final clustering labels
if (opt$a) {
  clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "allen" = allen,
                        "seurat" = seurat)
  } else {
  clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "seurat" = seurat)
}


# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
cat("Finished Consensus Merge\n")
saveRDS(object = mergers, file = paste0(output, "_mergers.rds"))