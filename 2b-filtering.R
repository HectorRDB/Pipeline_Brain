library(SingleCellExperiment)
library(stringr)

loc <- "/scratch/users/singlecell/MiniAtlas/data/"
dataset <- "SMARTer_nuclei_MOp/"
# source("1-loadData.R")

sce <- readRDS(file = paste0(loc, "rds/", str_replace(dataset, "/", ""), ".rds"))

sce <- lapply(sce, function(Sce){
  counts <- assays(Sce)$counts
  counts[is.na(counts)] <- 0
  assays(Sce) <- list(counts = counts, logcounts = logcounts(Sce))
  
  filt <- apply(assays(Sce)$counts, 1, function(x) {
    sum(x >= 50) >= 50
  })
  
  Sce <- Sce[filt, ]
  Sce
})

saveRDS(sce, file = paste0(loc, "rds/", str_replace(dataset, "/", ""), "_filt.rds"))
# in case fails in next steps...