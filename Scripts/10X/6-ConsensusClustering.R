library(optparse)

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1, type = "integer",
              help = "Number of cores to use"
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

library(clusterExperiment)
library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr)

# Load Data and clean seurat ----
# Load sc3 clustering results
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
k <- names(metadata(sc3)$sc3$consensus)
sc3 <- colData(sc3)[, paste0("sc3_", k, "_clusters")]
rm(k)

# Load RSEC and allen clustering results
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
allen <- colData(Rsec)$allenClusters
Rsec <- primaryCluster(Rsec)

# Load all seurat results and keep the two extrems
seurat <- readRDS(paste0(loc, "_seurat.rds"))
source("/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/scripts/10X/6-helper.R")
ARIs <- apply(seurat, 2, function(x) {
  apply(seurat, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

p <- plotARIs(ARIs, small = F) +
  ggtitle("Seurat concordance: ARIs for every pair of pais of parameters
          (resolution and k.param)")
ggsave(paste0(output, "_seurat_ARI.pdf"), p)

seurat_p <- c("0.6,50", "1.6,50")

seurat <- seurat[, seurat_p]
colnames(seurat) <- c("seurat1", "seurat2")

# Get the final clustering labels
clusMat <- data.frame("sc3" = sc3, "Rsec" = Rsec, "allen" = allen) %>%
  cbind(seurat)

# Do the consensus clustering ----
InitialARI <- apply(clusMat, 2, function(x) {
  apply(clusMat, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

p <- plotARIs(InitialARI) +
  ggtitle("ARI before any merging")
ggsave(paste0(output, "_Initial_ARI.pdf"), p)

clusMat2 <- clusMat[clusMat$Rsec != -1, ]
InitialARI2 <- apply(clusMat2, 2, function(x) {
  apply(clusMat2, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})
p <- plotARIs(InitialARI2) +
  ggtitle("ARI before any merging, no unclustered cells")
ggsave(paste0(output, "_Initial_ARI_no_unclus.pdf"), p)


print(paste0("Number of cores: ", opt$n))
mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)

FinalARI <- apply(mergers$currentMat, 2, function(x) {
  apply(mergers$currentMat, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

saveRDS(mergers, file = paste0(output, "_Consensus_Clustering.rds"))

p <- plotARIs(FinalARI) +
  ggtitle("ARI after merging")
ggsave(paste0(output, "_Final_ARI.pdf"), p)

clusMat2 <- mergers$currentMat[mergers$currentMat[,"Rsec"] != -1, ]

FinalARI2 <- apply(clusMat2, 2, function(x) {
  apply(clusMat2, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

p <- plotARIs(FinalARI2) +
  ggtitle("ARI after merging, no unclustered cells")
ggsave(paste0(output, "_Final_ARI_no_unclus.pdf"), p)
