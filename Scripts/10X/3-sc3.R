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

# Run sc3 per se ----
library(SC3)
library(scater)
library(stringr)
library(SingleCellExperiment)

# Add a normalization step ? 
sce <- readRDS(file = loc)
print("Dataset of size")
print(dim(sce))

rowData(sce)$feature_symbol <- rownames(sce)
# sce <- sc3_estimate_k(sce)
# K <- metadata(sce)$sc3$k_estimation
K <- c(60, 80, 100, 120, 140)
cat("Running the sc3 on a reduced set of ", round(.1 * ncol(sce)), "cells\n")
sce <- sc3(sce, ks = K, svm_max = ncol(sce) + 1, biology = FALSE,
           n_cores = as.numeric(opt$n), svm_num_cells = round(.1 * ncol(sce)))
cat("Fitting the other cells", "\n")
sce <- sc3_run_svm(sce, ks = K)

print(paste0("Saving output at ", output))
saveRDS(sce, file = output)
