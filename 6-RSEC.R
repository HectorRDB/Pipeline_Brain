library(clusterExperiment)
library(stringr)
library(zinbwave)
library(optparse)
library(SummarizedExperiment)

option_list <- list(
  make_option(c("-d", "--dataset"),
              action = "store", default = NA, type = "character",
              help = "The dataset on which to run the script"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$d)) {
  dataset <- opt$d
  cat("The selected dataset is", dataset)
} else {
  stop("A dataset needs to be provided")
}

# dataset <- "SMARTer_cells_MOp/"

loc <- "/scratch/users/singlecell/MiniAtlas/data/"
zinb <- readRDS(file = paste0(loc, "rds/", str_replace(dataset, "/", ""),
                             "_zinbWs.rds"))
sce <- readRDS(file = paste0(loc, "rds/", str_replace(dataset, "/", ""),
                             "_filt.rds"))

sequential <- FALSE
subsample <- T
clusterFunction <- "pam"

NCORES <- 8

# Add reduced dim
for (name in names(sce)) {
  Sce <- sce[[name]]
  for (i in 1:length(zinb[[name]])) {
    type <- names(zinb[[name]])[i]
    reducedDim(Sce, type = type) <- zinb[[name]][[i]]
  }
  sce[[name]] <- Sce
}
rm(Sce, type, i, name)

##---- run RSEC
for (name in names(sce)) {
  print(name)
  Sce <- sce[[name]]
  print(system.time(
    Sce <- RSEC(Sce, k0s = seq(5, 50, by = 5), alphas = c(0.1,0.3),
                reduceMethod = "zinb100",
                nReducedDims = c(3, 5, 10, 20, 50, 100),
                sequential = sequential, subsample = subsample, minSizes = 1,
                betas = c(0.8), clusterFunction = clusterFunction,
                ncores = NCORES, isCount = FALSE, dendroReduce = "zinb100",
                dendroNDims = 100, verbose = TRUE, consensusProportion = 0.7,
                subsampleArgs = list(resamp.num = 50, clusterFunction = "kmeans"),
                mergeMethod = "adjP", mergeCutoff = 0.1, mergeLogFCcutoff = 1,
                random.seed = 23578, consensusMinSize = 10, run = TRUE)
  ))
  sce[[name]] <- Sce
  rm(Sce)
}


saveRDS(sce, file = paste0(loc, "rds/", str_replace(dataset, "/", ""), "_RSEC.rds"))