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

# Filter data per se ----
 
# source("1-loadData.R")
sce <- readRDS(file = loc)

counts <- assays(sce)$counts
counts[is.na(counts)] <- 0
assays(sce) <- list(counts = counts, logcounts = logcounts(sce))

filt <- apply(assays(sce)$counts, 1, function(x) {
  sum(x >= 50) >= 50
})
  
sce <- sce[filt, ]

print(cat("Saving output at ", output))
saveRDS(sce, file = output)