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
library(SummarizedExperiment)

sce <- readRDS(file = loc)

clusterMatrixs <- list()
# Setup ----
sSeurat <- CreateSeuratObject(raw.data = assays(sce)$counts, project = 'allen40K')
sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
sSeurat <- FindVariableGenes(object = sSeurat, mean.function = ExpMean,
                             dispersion.function = LogVMR, do.plot = T)
sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nUMI")
sSeurat <- RunPCA(object = sSeurat, pc.genes = sSeurat@var.genes, 
                  do.print = TRUE, pcs.compute = 50, pcs.print = 1:5,
                  genes.print = 5)
  
# Run clustering ----
clusterMatrix <- NULL
for (RESOLUTION in c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)) {
  for (K.PARAM in c(30, 50, 100)) {
    sSeurat_star <- FindClusters(object = sSeurat, reduction.type = "pca",
                                 dims.use = 1:50, resolution = RESOLUTION,
                                 print.output = 0, k.param = K.PARAM,
                                 save.SNN = TRUE)
    clusterMatrix <- cbind(clusterMatrix, sSeurat_star@ident)
    colnames(clusterMatrix)[ncol(clusterMatrix)] <- paste(RESOLUTION, K.PARAM,
                                                          sep = ",")
  }
}

saveRDS(clusterMatrix, file = output)
