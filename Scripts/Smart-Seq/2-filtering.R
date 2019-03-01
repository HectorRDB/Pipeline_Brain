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
  make_option(c("-c", "--cell-cutoff"),
              action = "store", default = 50, type = "integer",
              help = "Which cutoff to use when filtering for cell"
  ),
  make_option(c("-r", "--read-cutoff"),
              action = "store", default = 50, type = "integer",
              help = "Which cutoff to use when filtering for read"
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
library(dplyr)
library(stringr)
library(SingleCellExperiment)

sce <- readRDS(file = loc)

counts <- assays(sce)$counts
counts[is.na(counts)] <- 0
assays(sce) <- list(counts = counts, logcounts = logcounts(sce))

print(opt$c)
print(opt$r)
filt <- rowSums(counts(sce) >= opt$r) >= opt$c
  
sce <- sce[filt, ]
dim(sce)
mean(counts(sce) == 0)

print(cat("Saving output at ", output))
saveRDS(sce, file = output)