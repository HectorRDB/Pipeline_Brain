suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-p", "--plot-output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the plots"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1, type = "integer",
              help = "Number of cores to use"
  ),
  make_option(c("-a", "--allen"),
              action = "store", default = T, type = "logical",
              help = "Wether to use allen or not"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc)
} else {
  stop("Missing l argument")
}
if (!is.na(opt$p)) {
  output_p <- opt$p
} else {
  stop("Missing p argument")
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
library(tidyverse)

# Load Data and clean seurat ----
# Load sc3 clustering results
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
k <- names(metadata(sc3)$sc3$consensus)
sc3 <- colData(sc3)[, paste0("sc3_", k, "_clusters")]
rm(k)

# Load RSEC and allen clustering results
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
allen <- colData(Rsec)$allenClusters
RsecT <- assignUnassigned(Rsec, clusterLabel = "Assigned")
RsecT <- primaryCluster(RsecT)
Rsec <- primaryCluster(Rsec)

# Load all seurat results and keep the two extrems
seurat <- readRDS(paste0(loc, "_seurat.rds"))
source("/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Scripts/Smart-Seq/7-helper.R")
# ARIs <- apply(seurat, 2, function(x) {
#  apply(seurat, 2, function(y) {
#     mclust::adjustedRandIndex(x, y)
#   })
# })

# p <- plotARIs(ARIs, small = F) +
#   ggtitle("Seurat concordance: ARIs for every pair of pais of parameters
#           (resolution and k.param)")
# ggsave(paste0(output_p, "_seurat_ARI.pdf"), p)

seurat_p <- "1.6,50"

seurat <- seurat[, seurat_p]

# Get the final clustering labels
if (opt$a) {
  clusMat <- data.frame("sc3" = sc3, "Rsec" = Rsec, "allen" = allen,
                        "seurat" = seurat)
  clusMatT <- data.frame("sc3" = sc3, "RsecT" = RsecT, "allen" = allen,
                         "seurat" = seurat)
  } else {
  clusMat <- data.frame("sc3" = sc3, "Rsec" = Rsec, "seurat" = seurat)
  clusMatT <- data.frame("sc3" = sc3, "RsecT" = RsecT, "seurat" = seurat)
}

# Inital plots ----
# No unclustered cells for RSEC
InitialARI <- apply(clusMat, 2, function(x) {
  apply(clusMat, 2, function(y) {
    inds <- x != -1 & y != -1
    xa <- x[inds]
    ya <- y[inds]
    mclust::adjustedRandIndex(xa, ya)
  })
})

p <- plotARIs(InitialARI) +
  ggtitle("ARI before any merging, no unclustered cells for RSEC")
ggsave(paste0(output_p, "_Initial_ARI.pdf"), p)

# All cells assigned
InitialARI <- apply(clusMatT, 2, function(x) {
  apply(clusMatT, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

p <- plotARIs(InitialARI) +
  ggtitle("ARI before any merging, all cells assigned")
ggsave(paste0(output_p, "_Initial_ARI_all_cells_Assigned.pdf"), p)


# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
cat("Finished Consensus Merge\n")
mergers$initalMat <- cbind(mergers$initalMat, clusMat[,"RsecT"])
colnames(mergers$initalMat)[4 + opt$a] <- "RsecT"
saveRDS(object = mergers, file = paste0(output, "_mergers.rds"))

# Final plots ----
# How is the number of cells reduced
pre <- apply(clusMat, 2, function(x) length(unique(x)))
post <- apply(mergers$currentMat, 2, function(x) length(unique(x)))
df <- data.frame(methods = names(pre),
                 before = pre,
                 after = post) %>%
  gather(key = "time", value = "Nb", -methods) %>%
  mutate(time = factor(time, levels = c("before", "after")))
p <- ggplot(df, aes(x = methods, y = Nb, fill = time)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_viridis_d(option = "E") +
  labs(x = "Clustering Methods", y = "Number of clusters", fill = "") +
  ggtitle("Reduction in number of clusters with consensus merging")
ggsave(paste0(output_p, "_clusters_reduction.pdf"), p)

# How is ARI improved
## Rsec with all cells
currentMat <- mergers$currentMat
Rsec_merges <- mergers$merges
Rsec_merges <- Rsec_merges[Rsec_merges[,1] == 2, ]
Rsec_merges <- Rsec_merges[, -1]
currentMat[, "Rsec"] <- lapply(1:nrow(currentMat), function(i) {
    cell <- clusMat[i ,"Rsec"]
    cellT <- clusMatT[i ,"RsecT"]
    if (cell == -1) {
      for (j in 1:nrow(Rsec_merges)) {
        if (cellT %in% Rsec_merges[j, ]) cellT <- min(Rsec_merges[j, ])
      }
    }
    return(cellT)
  }) %>%
  unlist()

FinalARI <- apply(currentMat, 2, function(x) {
  apply(currentMat, 2, function(y) {
    mclust::adjustedRandIndex(x, y)
  })
})

p <- plotARIs(FinalARI) +
  ggtitle("ARI after merging, all cells assigned")
ggsave(paste0(output_p, "_Final_ARI_all_cells_assigned.pdf"), p)

## No unclustered cells for Rsec
FinalARI <- apply(mergers$currentMat, 2, function(x) {
  apply(mergers$currentMat, 2, function(y) {
    inds <- x != -1 & y != -1
    xa <- x[inds]
    ya <- y[inds]
    mclust::adjustedRandIndex(xa, ya)
  })
})

p <- plotARIs(FinalARI) +
  ggtitle("ARI after merging, no unclustered cells for RSEC")
ggsave(paste0(output_p, "_Final_ARI.pdf"), p)
