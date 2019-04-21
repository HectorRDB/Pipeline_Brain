suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-p", "--plot"),
              action = "store", default = NA, type = "character",
              help = "Where to store the graphical output"
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

library(clusterExperiment)
library(stringr)
library(zinbwave)
library(SingleCellExperiment)

# Load data ----
sce <- readRDS(file = loc)

##---- run RSEC
sequential <- FALSE
subsample <- TRUE
clusterFunction <- "pam"
NCORES <- as.numeric(opt$n)
reduceMeth <- reducedDimNames(sce)

cat("Cluster Many\n")

print(system.time(
  sce <- clusterMany(sce,
                     ks = seq(10, 50, by = 5),
                     clusterFunction = clusterFunction,
                     reduceMethod = reduceMeth,
                     alphas = c(0.1, 0.3),
                     sequential = sequential,
                     betas = c(0.8),
                     minSizes = 1,
                     isCount = TRUE,
                     subsampleArgs = list(resamp.num = 50,
                                          clusterFunction = "kmeans"),
                     random.seed = 23578,
                     run = TRUE,
                     ncores = NCORES)
))
saveRDS(sce, file = output)

cat("make Consensus\n")
print(system.time(
  sce <- makeConsensus(sce,
                       minSize = 10,
                       proportion = 0.7)
))
saveRDS(sce, file = output)

cat("make Dendogram\n")
print(system.time(
  sce <- makeDendrogram(sce,
                        reduceMethod = reduceMeth[length(reduceMeth)],
                        nDims = 50)
))
saveRDS(sce, file = output)

cat("Merge Clusters\n")
print(system.time(
  sce <- mergeClusters(sce,
                       mergeMethod = "adjP",
                       cutoff = 0.05,
                       logFCcutoff = 1)
))

saveRDS(sce, file = output)

pdf(file = opt$p)
plotClusters(sce)
dev.off()

