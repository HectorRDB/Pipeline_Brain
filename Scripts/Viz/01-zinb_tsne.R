suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the object after running"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-K", "--Kvalue"),
              action = "store", default = 10,
              help = "The value of K [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc, "\n")
} else {
  stop("Missing l argument\n")
}
if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument\n")
}
if (!is.na(opt$K)) {
  K <- opt$K
} else {
  stop("Missing K argument\n")
}

# Loading the data ----
library(stringr)
library(vctrs)
library(dplyr)
library(BiocParallel)
library(Rtsne)
library(SingleCellExperiment)
library(zinbwave)
library(matrixStats)
library(ggplot2)

if (str_detect(loc, "SMARTer")) {
  sce <- readRDS(loc)
  zinbW <- reducedDims(sce)[[paste0("zinb-K", K)]]
  rownames(zinbW) <- colnames(sce)
} else {
  zinbW <- readRDS(loc)
}
TNSE <- Rtsne(zinbW, initial_dims = min(50, K))

write.csv(data.frame(cells = rownames(zinbW), x = TNSE$Y[, 1], y = TNSE$Y[, 2]),
          file = output)
