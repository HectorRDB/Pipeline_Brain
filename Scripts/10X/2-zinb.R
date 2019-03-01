suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
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
  ),
  make_option(c("-i", "--initial"),
              action = "store", default = 10,
              help = "The smallest value of K [default %default]"
  ),
  make_option(c("-f", "--final"),
              action = "store", default = 50,
              help = "The largest value of K [default %default]"
  ),
  make_option(c("-d", "--dims"),
              action = "store", default = 9,
              help = "How many values of K between i and f [default %default]"
  ),
  make_option(c("-t", "--tests"),
              action = "store", default = F,
              help = "Make plots to check validity [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc, "\n")
} else {
  stop("Missing l argument\n")
}
if (!is.na(opt$r)) {
  output_r <- opt$r
} else {
  stop("Missing r argument\n")
}

library(stringr)
library(clusterExperiment)
library(dplyr)
library(BiocParallel)
library(zinbwave)
library(matrixStats)

# Load data ----
sce <- readRDS(file = loc)

# Run ZinbWave ----
NCORES <- as.numeric(opt$n)
BiocParallel::register(MulticoreParam(NCORES))

vars <- matrixStats::rowVars(logcounts(sce))
cat("Running with K = 0 on the full data\n")
cat("Number of cores:", NCORES, "\n")
cat("Time to run zinbwave (seconds):\n")
print(system.time(zinb0 <- zinbwave(sce)))

ind <- vars > sort(vars,decreasing = TRUE)[1000]
whichGenes <- rownames(sce)[ind]
zinbDims <- floor(seq(from = opt$i, to = opt$f, length.out = opt$d))
cat("Using the following values for K :", zinbDims)
sceVar <- sce[ind,]

zinbWs <- lapply(zinbDims, function(zinbDim) {
  cat("Running with K = ", zinbDim, " on the filtered data\n")
  cat("Number of cores:", NCORES, "\n")
  cat("Time to run zinbwave (seconds):\n")
  print(system.time(zinb <- zinbwave(sceVar, K = zinbDim)))
  return(zinb)
})

for (i in 1:length(zinbWs)) {
  type <- paste0("zinb-K", zinbDims[i])
  reducedDim(sce, type = type) <- reducedDim(zinbWs[[i]])
}
saveRDS(sce, file = paste0(output_r, ".rds"))

if (opt$t) {
  library(ggplot2)
  library(Rtsne)
  
  clusters <- read.csv(
    "/scratch/users/singlecell/MiniAtlas/data/SMARTer_nuclei_MOp/cluster.annotation.csv",
    header = T)
  cols2 <- clusters$cluster_color %>% as.character()
  names(cols2) <- clusters$cluster_id %>% as.character()
  
  for (i in 1:length(zinbWs)) {
    zinbW <- reducedDim(zinbWs[[i]])
    cat("Tsne plot for dim", zinbDims[i])
    TNSE <- Rtsne(zinbW)
    print("....tsne")
    df <- data.frame(x = TNSE$Y[, 1], y = TNSE$Y[, 2],
                     cols = as.factor(colData(sce)$allenClusters))
    p <- ggplot(df, aes(x = x, y = y, col = cols)) +
      geom_point() +
      theme_classic() +
      scale_color_manual(values = cols2, breaks = names(cols2)) +
      labs(x = "dim1", y = "dim2")
    ggsave(paste0(output_r, "_K", zinbDims[i], ".pdf"),p)
    cat("....saving plot")
  }
}
