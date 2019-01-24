library(Seurat)
library(dplyr)
library(clusterExperiment)

dataset <- "SMARTer_cells_MOp/"
source("3a-zinb.R")

# Setup ----
sSeurat <- CreateSeuratObject(raw.data = assays(sce)$counts, project = 'allen40K')

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

saveRDS(clusterMatrix, file = 'seurat_manyclusters.rds')
