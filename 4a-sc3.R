library(SC3)
library(scater)
library(stringr)

dataset <- "SMARTer_cells_MOp/"
loc <- "/scratch/users/singlecell/MiniAtlas/data/"
# source("2-filtering.R")
sce <- readRDS(file = paste0(loc, "rds/", str_replace(dataset, "/", ""),
                             "_filt.rds"))
NCORES <- 8

for (name in names(sce)) {
  Sce <- sce[[name]]
  rowData(Sce)$feature_symbol <- rownames(Sce)
  Sce <- sc3_estimate_k(Sce)
  K <- metadata(Sce)$sc3$k_estimation
  save(K, file = paste0(loc, "RData/", str_replace(dataset, "/", ""),
                        "_sc3_k.RData"))
  Sce <- sc3(Sce, ks = K, svm_max = ncol(Sce) + 1, biology = FALSE,
             n_cores = NCORES)
  sce[[name]] <- Sce
  print(name)
}

saveRDS(sce, file = paste0(loc, "rds/", str_replace(dataset, "/", ""), "sc3_out.rds"))


