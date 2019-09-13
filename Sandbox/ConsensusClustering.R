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

# Load Data----
library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Load Monocle clustering results
# Load sc3  and allen clustering results
print("Loading sc3")
sc3 <- readRDS(paste0(loc, "_sc3.rds"))
sc3 <- data.frame("sc3_100_clusters" = colData(sc3)[, "sc3_100_clusters"],
                  cells = colnames(sc3))
write.csv(sc3, file = output)