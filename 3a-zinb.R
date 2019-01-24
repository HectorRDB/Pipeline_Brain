# dataset <- "SMARTer_cells_MOp/"
source("2-filtering.R")

library(clusterExperiment)
library(dplyr)
library(BiocParallel)
library(zinbwave)

# Subset ---- 
sce_small <- sce
# nvars <- 5000
# vars <- rowVars(assays(sce)$logcounts)
# names(vars) <- rownames(sce)
# vars <- sort(vars, decreasing = TRUE)
# sce_small <- sce[names(vars)[1:nvars],] #only the most variable genes.
# sce_small <- sce_small[, colSums(assays(sce)$counts) > 0]

# Run ZinbWave ----
library(BiocParallel)
library(zinbwave)
NCORES <- 8
BiocParallel::register(MulticoreParam(NCORES))

zinbWs <- list()
for (zinbDim in c(3, 5, 10, 20, 50, 100)) {
  cat("Number of cores:", NCORES, "\n")
  cat("Time to run zinbwave (seconds):\n")
  print(system.time(zinb <- zinbwave(sce_small, K = zinbDim)))
  zinbWs[[length(zinbWs) + 1]] <- reducedDim(zinb)
}
names(zinbWs) <- paste0("zinb", c(3, 5, 10, 20, 50, 100))

# saveRDS(zinbWs, file = "zinbWs.rds") 
# in case fails in next steps...

