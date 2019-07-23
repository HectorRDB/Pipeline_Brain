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
sc3 <- colData(sc3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
rm(k)

# Load monocle and allen clustering results
Monocle <- readRDS(paste0(loc, "_monocle.rds"))
Names <- colnames(Monocle)
allen <- pData(Monocle)$allenClusters %>% as.numeric()
Monocle <- pData(Monocle)$Cluster %>% as.numeric()

# Load RSEC results
Rsec <- readRDS(paste0(loc, "_RSEC.rds"))
Rsec <- mergeClusters(Rsec,
                      mergeMethod = "adjP",
                      plotInfo = "adjP",
                      cutoff = 0.01,
                      clusterLabel = "Clusters",
                      plot = F,
                      DEMethod = "limma")

# Load all seurat results and keep one of them
seurat <- readRDS(paste0(loc, "_seurat.rds"))
source("/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Scripts/Smart-Seq/8-helper.R")

seurat_p <- "1.6,50"

seurat <- seurat[, seurat_p] %>% as.numeric()
