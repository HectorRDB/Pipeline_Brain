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

suppressMessages(library(reticulate))
suppressMessages(library(monocle))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
suppressMessages(library(zinbwave))
import("louvain")

# Load data and convert to Delayed Array ----
sce <- readRDS(file = paste0(loc, "_monocle.rds"))
zinbW <- readRDS(file = paste0(loc, "_zinb.rds"))
sce <- newCellDataSet(sce@assayData$exprs,
                      phenoData = sce@phenoData,
                      featureData = sce@featureData)

# Pre-process
sce@normalized_data_projection <- zinbW
sce@auxOrderingData$normalize_expr_data <- sce@assayData$exprs
print("Doing the reduced dimension")
sce <- reduceDimension(sce,
                       max_components = 2,
                       reduction_method = 'UMAP',
                       metric = "correlation",
                       min_dist = 0.75,
                       n_neighbors = 50,
                       verbose = T)

# run Monocle ----
print("Running Monocle")
print(system.time(
  sce <- clusterCells(sce,
                      method = 'louvain',
                      res = 1e-2,
                      louvain_iter = 1,
                      verbose = T)
))

saveRDS(sce, file = output)
