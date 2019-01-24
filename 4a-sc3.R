library(clusterExperiment)
library(SC3)
library(scater)

library(clusterExperiment)

dataset <- "SMARTer_cells_MOp/"
source("3a-zinb.R")

sce <- SingleCellExperiment(assays = list(counts = assays(sce)$counts,
                                          logcounts = assays(ce)$logcounts),
                            colData = colData(sce), rowData = rowData(sce))
rowData(sce)$feature_symbol <- rownames(sce)

sce <- sc3_estimate_k(sce)
K <- metadata(sce)$sc3$k_estimation
save(K, file = "inhibit/sc3_k.RData")
NCORES <- 8

sce <- sc3(sce, ks = K, svm_max = ncol(sce) + 1, biology = FALSE, n_cores = N)

saveRDS(sce, file = "inhibit/sc3_out.rds")

