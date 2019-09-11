suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output matrix"
  ),
  make_option(c("-s", "--sce"),
              action = "store", default = NA, type = "character",
              help = "Where to store the clusterExperiment Object"
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

if (!is.na(opt$s)) {
  clus <- opt$s
} else {
  stop("Missing s argument")
}

library(clusterExperiment)
library(stringr)
library(zinbwave)
library(SingleCellExperiment)

# Load data ----
sce <- readRDS(file = loc)

##---- run RSEC
sequential <- FALSE
subsample <- T
clusterFunction <- "pam"
NCORES <- as.numeric(opt$n)
reduceMeth <- reducedDimNames(sce)

print(system.time(
  sce <- RSEC(sce, k0s = seq(10, 50, by = 5), alphas = c(0.1, 0.3),
              reduceMethod = reduceMeth, sequential = sequential,
              subsample = subsample, minSizes = 1, betas = c(0.8), 
              clusterFunction = clusterFunction, ncores = NCORES, run = TRUE,
              isCount = FALSE, dendroReduce = reduceMeth[length(reduceMeth)],
              dendroNDims = 50, consensusProportion = 0.7, verbose = TRUE,
              random.seed = 23578, mergeMethod = "adjP", mergeCutoff = 0.05,
              subsampleArgs = list(resamp.num = 50, clusterFunction = "kmeans"),
              mergeLogFCcutoff = 1, consensusMinSize = 10)
))

# Saving objects
saveRDS(sce, clus)

RsecT <- assignUnassigned(sce, clusterLabel = "Assigned")
RsecT <- primaryCluster(RsecT) %>% as.numeric()
Rsec <- primaryCluster(sce) %>% as.numeric()
write.csv(data.frame("Rsec" = Rsec, "RsecT" = RsecT, cells = colnames(sce)),
          file = output)
