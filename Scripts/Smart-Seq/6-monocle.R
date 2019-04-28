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
sce <- readRDS(file = loc)
pd <- new("AnnotatedDataFrame", data = as.data.frame(sce@colData))
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(assays(sce)$counts)))
zinbW <- reducedDim(sce, type = reducedDimNames(sce)[3])
rownames(fd) <- rownames(assays(sce)$counts)
sce <- newCellDataSet(assays(sce)$counts,
                      phenoData = pd,
                      featureData = fd)

# Pre-process
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1005)
sce <- estimateSizeFactors(sce)
sce <- estimateDispersions(sce)
sce <- preprocessCDS(sce,
                     method = 'PCA',
                     norm_method = 'log',
                     num_dim = 50,
                     verbose = T)
sce <- reduceDimension(sce,
                       max_components = 2,
                       reduction_method = 'UMAP',
                       metric = "correlation",
                       min_dist = 0.75,
                       n_neighbors = 50,
                       verbose = T)

# run RSEC ----
print("Running RSEC")
print(system.time(
  sce <- clusterCells(sce,
                      method = 'louvain',
                      res = 1e-2,
                      louvain_iter = 1,
                      verbose = T)
))

# pData(sce)$Cluster
saveRDS(sce, file = output)
