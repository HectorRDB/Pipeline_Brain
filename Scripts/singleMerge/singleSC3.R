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
# .libPaths("/accounts/projects/epurdom/singlecell/R/x86_64-pc-linux-gnu-library/3.5")
print(.libPaths())
library(SC3)
library(scater)
library(stringr)
library(SingleCellExperiment)
library(purrr)

# Add a normalization step ? 
sce <- readRDS(file = loc)

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3_estimate_k(sce)
K <- metadata(sce)$sc3$k_estimation
ks <- 2:K
names(ks) <- ks
sc3 <- map_df(ks, function(k){
  SC3 <- sc3(sce, ks = k, svm_max = ncol(sce) + 1, biology = FALSE,
             n_cores = as.numeric(opt$n))
  SC3 <- colData(SC3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
  return(SC3)
})

sc3$cells <- colnames(sce)
write.table(sc3, file = output)
