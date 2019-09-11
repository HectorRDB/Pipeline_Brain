# Options for the script ----
suppressWarnings(library(optparse))
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the clustering files"
  ),
  make_option(c("-r", "--rsec"),
              action = "store", default = NA, type = "character",
              help = "Location of the Rsec object"
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
library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(tidyverse)
library(Dune)
library(mclust)
library(clusterExperiment)

# Load Data ----
# Load sc3 clustering results
sc3 <- read.csv(paste0(loc, "_SC3.csv"))
ggsave(filename = paste0(opt$p, "_monocle_ARI.png"),
       plot = clusterMatToAri(sc3 %>% select(-cells)))
Names <- sc3$cells
sc3 <- sc3[,"sc3_0_clusters"]

# Load Seurat clustering results
seurat <- read.csv(paste0(loc, "_seurat.csv"))
ggsave(filename = paste0(opt$p, "_seurat_ARI.png"),
       plot = clusterMatToAri(seurat %>% select(-cells)))
seurat_p <- "1.2,50"
seurat <- seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))
ggsave(filename = paste0(opt$p, "_monocle_ARI.png"),
       plot = clusterMatToAri(Monocle %>% select(-cells)))
monocle_p <- "k_45"
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Load RSEC clustering results
Rsec <- readRDS(opt$r)


for (clustering in c("sc3", "Monocle", "seurat")) {
  Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
}

# Doing the merges
cutoffs <- seq(from = 0, to = 1, by = .05)
res <- list()
for (clustering in c("sc3", "Monocle", "seurat", "Rsec")) {
  print(clustering)
  Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
  names(cutoffs) <- paste(clustering, cutoffs, sep = "_")
  res[[clustering]] <- map_df(cutoffs,
                function(i){
                  print(paste0("...", i))
                  Rsec3 <- mergeClusters(Rsec2,
                                         mergeMethod = "adjP",
                                         plotInfo = "adjP",
                                         cutoff = i,
                                         clusterLabel = "Clusters",
                                         plot = F,
                                         DEMethod = "limma")
                  return(Rsec3@clusterMatrix[,"Clusters"])  
                })
}

res <- do.call('cbind', res) %>% as.data.frame()
res$cells <- colnames(Rsec)
write_csv(res, path = output)