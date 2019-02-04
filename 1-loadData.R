library(dplyr)
library(stringr)
library(SingleCellExperiment)
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



# Load data per se ----

counts <- read.csv(paste0(loc, "exon.counts.csv.gz"),
                    header = T, row.names = 1)
colnames(counts) <- str_replace_all(colnames(counts), "\\.", "-")
meta <- read.csv(paste0(loc, "sample_metadata.csv.gz"),
                 header = T, row.names = 1)
allenClusters <- read.csv(paste0(loc, "cluster.membership.csv"),
                          header = T, col.names = c("sample", "clusters"))
meta$allenClusters <- allenClusters$clusters

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts),
                                          logcounts = as.matrix(log1p(counts))),
                            colData = meta)

cat("Saving output at ", output)
saveRDS(sce, file = output)
# in case fails in next steps...