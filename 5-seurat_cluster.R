library(optparse)
library(Seurat)
library(dplyr)
library(stringr)
library(SummarizedExperiment)

option_list <- list(
  make_option(c("-d", "--dataset"),
              action = "store", default = NA, type = "character",
              help = "The dataset on which to run the script"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$d)) {
  dataset <- opt$d
  cat("The selected dataset is", dataset)
} else {
  stop("A dataset needs to be provided")
}

# dataset <- "SMARTer_cells_MOp/"

loc <- "/scratch/users/singlecell/MiniAtlas/data/"
sce <- readRDS(file = paste0(loc, "rds/", str_replace(dataset, "/", ""),
                             "_filt.rds"))


clusterMatrixs <- list()
for (name in names(sce)) {
  Sce <- sce[[name]]
  # Setup ----
  sSeurat <- CreateSeuratObject(raw.data = assays(Sce)$counts, project = 'allen40K')
  
  sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
  
  ssSeurat <- FindVariableGenes(object = sSeurat, mean.function = ExpMean,
                                dispersion.function = LogVMR,
                                x.low.cutoff = 0.0125, x.high.cutoff = 3,
                                y.cutoff = 0.5)
  
  sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nUMI")
  
  sSeurat <- RunPCA(object = sSeurat, pc.genes = sSeurat@var.genes,
                    do.print = TRUE, pcs.compute = 50, pcs.print = 1:5,
                    genes.print = 5)
  
  # Run clustering ----
  clusterMatrix <- NULL
  for (RESOLUTION in c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)) {
    for (K.PARAM in c(10, 30, 50, 100)) {
      sSeurat_star <- FindClusters(object = sSeurat, reduction.type = "pca",
                                   dims.use = 1:50, resolution = RESOLUTION,
                                   print.output = 0, k.param = K.PARAM,
                                   save.SNN = TRUE)
      clusterMatrix <- cbind(clusterMatrix, sSeurat_star@ident)
      colnames(clusterMatrix)[ncol(clusterMatrix)] <- paste(RESOLUTION,
                                                            K.PARAM, sep = ",")
    }
  }
  clusterMatrixs[[name]] <- clusterMatrix
}

saveRDS(clusterMatrixs, file = paste0(loc, "rds/", str_replace(dataset, "/", ""),
                                      "_seurat_out.rds"))
