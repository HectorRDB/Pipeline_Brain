suppressWarnings(library(optparse))

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
  ),
  make_option(c("-a", "--allen"),
              action = "store", default = T, type = "logical",
              help = "Wether to use allen or not"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at ", loc, "\n")
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
  cat("The output will be stored at ", output, "\n")
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

# Load monocle and allen clustering results
Monocle <- readRDS(paste0(loc, "_monocle.rds"))
Names <- colnames(Monocle)
allen <- pData(Monocle)$allenClusters
Monocle <- pData(Monocle)$Cluster

# Load RSEC results
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
RsecT <- assignUnassigned(Rsec, clusterLabel = "Assigned")
RsecT <- primaryCluster(RsecT)
Rsec <- primaryCluster(Rsec)

# Load all seurat results and keep one of them
seurat <- readRDS(paste0(loc, "_seurat.rds"))
source("/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Scripts/Smart-Seq/8-helper.R")

seurat_p <- "1.6,50"

seurat <- seurat[, seurat_p]

# Get the final clustering labels
if (opt$a) {
  clusMat <- data.frame("sc3" = sc3, "Rsec" = Rsec, "Monocle" = Monocle,
                        "allen" = allen, "seurat" = seurat)
  clusMatT <- data.frame("sc3" = sc3, "RsecT" = RsecT, "Monocle" = Monocle,
                         "allen" = allen, "seurat" = seurat)
  } else {
  clusMat <- data.frame("sc3" = sc3, "Rsec" = Rsec, "Monocle" = Monocle,
                        "seurat" = seurat)
  clusMatT <- data.frame("sc3" = sc3, "RsecT" = RsecT, "Monocle" = Monocle,
                         "seurat" = seurat)
}
rownames(clusMat) <- Names  

# Do the consensus clustering ----
print(paste0("Number of cores: ", opt$n))
print(system.time(
  mergers <- mergeManyPairwise(clusteringMatrix = clusMat, nCores = opt$n)
))
cat("Finished Consensus Merge\n")
mergers$initalMat <- cbind(mergers$initalMat, clusMatT[,"RsecT"])
colnames(mergers$initalMat)[5 + opt$a] <- "RsecT"
saveRDS(object = mergers, file = paste0(output, "_mergers.rds"))