# Options for the script ----
suppressWarnings(library(optparse))
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
# Loading files ---- 
library(vctrs)
library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
.libPaths("/accounts/projects/epurdom/singlecell/R/x86_64-pc-linux-gnu-library/3.5")
library(clusterExperiment)
library(monocle)

# Load sc3 clustering results
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
k <- names(metadata(sc3)$sc3$consensus)
sc3 <- colData(sc3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
rm(k)

# Load monocle clustering results
Monocle <- readRDS(paste0(loc, "_monocle.rds"))
monocle_p <- "k_45"
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Load RSEC results
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
Rsec <- assignUnassigned(Rsec, clusterLabel = "Rsec")

# Load all seurat results and keep one of them
seurat <- readRDS(paste0(loc, "_seurat.rds"))
seurat_p <- "1.2,50"
seurat <- seurat[, seurat_p] %>% as.numeric()

for (clustering in c("sc3", "Monocle", "seurat")) {
  Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
}

# Doing the merges
cutoffs <- seq(from = .05, to = 1, by = .05)
res <- list()
for (clustering in c("sc3", "Monocle", "seurat", "Rsec")) {
  print(clustering)
  Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
  names(cutoffs) <- paste(clustering, cutoffs, sep = "_")
  res[[clustering]] <- map_df(cutoffs,
                function(i){
                  print(paste0("...", i))
                  Rsec2 <- mergeClusters(Rsec2,
                                         mergeMethod = "adjP",
                                         plotInfo = "adjP",
                                         cutoff = i,
                                         clusterLabel = "Clusters",
                                         plot = F,
                                         DEMethod = "limma")  
                  return(Rsec2@clusterMatrix[,"Clusters"])  
                })
}

res <- do.call('cbind', res)
write_csv(res, path = output)