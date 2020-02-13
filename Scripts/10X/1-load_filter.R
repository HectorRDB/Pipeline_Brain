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
library(Matrix)


# Load data per se ----

cat("Loading the data", "\n")
if (str_detect(loc, "MOp")) {
  print("This is an Allen dataset")
  counts <- Read10X_h5(paste0(loc, "umi_counts.h5"))
  
  colnames(counts) <- str_replace_all(colnames(counts), "\\.", "-")
  print("Getting the metadata")
  meta <- read.csv(paste0(loc, "sample_metadata.csv"),
                   header = T, row.names = 1)
  allenClusters <- read.csv(paste0(loc, "cluster.membership.csv"),
                            header = T, col.names = c("sample", "clusters"),
                            stringsAsFactors = FALSE)
  meta <- meta[allenClusters$sample, ]
  meta$allenClusters <- allenClusters$clusters
  meta <- meta %>%
    mutate(sample = str_extract(exp_component_name, "^[A-Z]*-[0-9]*"))
} else {
  print("This is a dataset from someone else, assuming a csv input file")
  counts <- read.csv(loc, row.names = 1)
  meta <- data.frame(cells = colnames(counts))
}

print("Preparing the data")
print("Initially, we have ", ncol(counts), " samples and ", nrow(counts),
      "genes.")
counts <- counts[, meta$sample]
print("Then, we keep ", ncol(counts), " good quality cells")
print(quantile(counts, p = (0:10)/10))
filt <- rowSums(counts >= opt$c) >= opt$c
print(sum(filt))
print(sum(counts[filt, ]) / sum(counts))
counts <- counts[filt, ]
print("At the end, we have ", ncol(counts), " samples and ", nrow(counts),
      "genes.")

cat("Saving output to ", output)
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts),
                                          logcounts = as.matrix(log1p(counts))),
                            colData = meta)

saveRDS(sce, file = output)
