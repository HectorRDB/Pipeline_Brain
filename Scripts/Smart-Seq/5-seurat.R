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

# Load data ----
library(Seurat)
library(dplyr)
library(stringr)
library(SingleCellExperiment)

sce <- readRDS(file = loc)

# Setup ----
sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K')
sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                dispersion.function = LogVMR, do.plot = T)
sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nCount_RNA")
sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 100)

# Run clustering ----
clusterMatrix <- NULL
for (RESOLUTION in seq(from = 0.3, to = 2.5, by = .1)) {
  print(RESOLUTION)
  for (K.PARAM in c(30, 50, 100)) {
    print(paste0("...", K.PARAM))
    sSeurat_star <- FindNeighbors(sSeurat, dims = 1:K.PARAM)
    sSeurat_star <- FindClusters(sSeurat_star, resolution = RESOLUTION)
    clusterMatrix <- cbind(clusterMatrix, Idents(sSeurat_star))
    colnames(clusterMatrix)[ncol(clusterMatrix)] <- paste(RESOLUTION, K.PARAM,
                                                          sep = ",")
  }
}

clusterMatrix <- as.data.frame(clusterMatrix)
clusterMatrix$cells <- colnames(sce)
write.csv(clusterMatrix, file = output)
