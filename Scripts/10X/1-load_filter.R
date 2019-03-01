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
  make_option(c("-c", "--cutoff"),
              action = "store", default = 50, type = "character",
              help = "The cutoff for filtering"
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

library(dplyr)
library(stringr)
library(Seurat)
library(SingleCellExperiment)

# Load data per se ----

cat("Loading the data", "\n")
counts <- Read10X_h5(paste0(loc, "umi_counts.h5"))

colnames(counts) <- str_replace_all(colnames(counts), "\\.", "-")
meta <- read.csv(paste0(loc, "sample_metadata.csv"),
                 header = T, row.names = 1)
allenClusters <- read.csv(paste0(loc, "cluster.membership.csv"),
                          header = T, col.names = c("sample", "clusters"))
meta$allenClusters <- allenClusters$clusters

cat("Preparing the data", "\n")


counts[is.na(counts)] <- 0
counts <- as.matrix(counts)
filt <- rowSums(counts >= opt$c) >= opt$c
cat(sum(counts[filt, ]) / sum (counts))
counts <- counts[filt, ]

cat(mean(filt), "\n")
cat(sum(filt), "\n")

cat("Saving output to ", output)

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts),
                                          logcounts = as.matrix(log1p(counts))),
                            colData = meta)

saveRDS(sce, file = output)
# # in case fails in next steps...
