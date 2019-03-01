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

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3_estimate_k(sce)
K <- metadata(sce)$sc3$k_estimation
sce <- sc3(sce, ks = K, svm_max = ncol(sce) + 1, biology = FALSE,
           n_cores = as.numeric(opt$n))

print(cat("Saving output at ", output))
saveRDS(sce, file = output)
