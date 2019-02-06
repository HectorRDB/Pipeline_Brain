library(optparse)

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output-normalized"),
              action = "store", default = NA, type = "character",
              help = "Where to store the normalized output"
  ),
  make_option(c("-r", "--reduced-dim-output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the reduced dim object"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1,
              help = "Number of cores to use [default %default]"
  )
)

library(stringr)
library(clusterExperiment)
library(dplyr)
library(BiocParallel)
library(zinbwave)
library(matrixStats)

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

# Load data ----
sce <- readRDS(file = loc)

# Run ZinbWave ----
NCORES <- as.numeric(opt$n)
BiocParallel::register(MulticoreParam(NCORES))

vars <- matrixStats::rowVars(logcounts(sce))
ind <- vars > sort(vars,decreasing = TRUE)[1000]
whichGenes <- rownames(sce)[ind]

zinbDims <- 10 * 1:5
cat("Running with K = 0 on the full data\n")
cat("Number of cores:", NCORES, "\n")
cat("Time to run zinbwave (seconds):\n")
print(system.time(zinb0 <- zinbwave(sce)))

sceVar <- sce[ind,]

success <- lapply(zinbDims, function(zinbDim) {
  cat("Running with K = ", zinbDim, " on the filtered data\n")
  cat("Number of cores:", NCORES, "\n")
  cat("Time to run zinbwave (seconds):\n")
  print(system.time(zinb <- zinbwave(sceVar, K = zinbDim)))
  type <- zinbDim
  reducedDim(sceVar, type = paste0("zinb", type)) <- zinb
  cat("Saving output at ", output, "\n")
  save(c(sceVar, zinb0), file = output)
  return(zinbDim)
})

cat("Saving output at ", output, "\n")
save(c(sceVar, zinb0), file = output)