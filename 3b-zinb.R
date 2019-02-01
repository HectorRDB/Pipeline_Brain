dataset <- "SMARTer_nuclei_MOp/"
loc <- "/scratch/users/singlecell/MiniAtlas/data/"
library(stringr)
library(clusterExperiment)
library(dplyr)
library(BiocParallel)
library(zinbwave)

# source("2-filtering.R")
sce <- readRDS(file = paste0(loc, "rds/", str_replace(dataset, "/", ""),
                             "_filt.rds"))

# Subset ---- 

# sce_small <- sce
# nvars <- 5000
# vars <- rowVars(assays(sce)$logcounts)
# names(vars) <- rownames(sce)
# vars <- sort(vars, decreasing = TRUE)
# sce_small <- sce[names(vars)[1:nvars],] #only the most variable genes.
# sce_small <- sce_small[, colSums(assays(sce)$counts) > 0]

# Run ZinbWave ----
NCORES <- 8
BiocParallel::register(MulticoreParam(NCORES))

zinbWs <- list(Inhibit = list(),
               Excite = list())
for (name in names(sce)) {
  sce_small <- sce[[name]]
  for (zinbDim in c(3, 5, 10, 20, 50, 100)) {
    cat("Number of cores:", NCORES, "\n")
    cat("Time to run zinbwave (seconds):\n")
    print(system.time(zinb <- zinbwave(sce_small, K = zinbDim)))
    zinbWs[[name]][[length(zinbWs[[name]]) + 1]] <- reducedDim(zinb)
  } 
  names(zinbWs[[name]]) <- paste0("zinb", c(3, 5, 10, 20, 50, 100))
}

saveRDS(zinbWs, file = paste0(loc, "rds/", str_replace(dataset, "/", ""), "_zinbWs.rds"))

# in case fails in next steps...

