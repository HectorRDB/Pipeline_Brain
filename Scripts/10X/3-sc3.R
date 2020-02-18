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
library(purrr)
library(SingleCellExperiment)

# Add a normalization step ? 
sce <- readRDS(file = loc)
rowData(sce)$feature_symbol <- rownames(sce)
# ks <- c(80, 100, 120)
ks <- 100
names(ks) <- ks
cat("Running the sc3 on a reduced set of ", round(.1 * ncol(sce)), "cells\n")

sc3 <- map_df(ks, function(k){
  SC3 <- sc3(sce, ks = k, svm_max = ncol(sce) + 1, biology = FALSE,
             n_cores = as.numeric(opt$n), svm_num_cells = round(.1 * ncol(sce)))
  SC3 <- sc3_run_svm(SC3, ks = k)
  SC3 <- colData(SC3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
  return(SC3)
})

sc3$cells <- colnames(sce)
write.csv(sc3, file = output)
