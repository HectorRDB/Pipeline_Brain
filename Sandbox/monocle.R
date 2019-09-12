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

suppressMessages(library(reticulate))
suppressMessages(library(monocle3))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(flexclust))
suppressMessages(library(purrr))
suppressMessages(library(mcclust))
suppressMessages(library(zinbwave))
import("louvain")

# Load data and convert----
sce <- readRDS(file = loc)
pd <- data.frame(cells = rownames(sce))
fd <- data.frame(gene_short_name = colnames(sce))
zinbW <- sce[,3]
rownames(fd) <- fd$gene_short_name
sce <- new_cell_data_set(t(sce),
                         cell_metadata = pd,
                         gene_metadata = fd)

# Pre-process
sce@reducedDims <- SimpleList("PCA" = zinbW)

print("Doing the reduced dimension")
sce <- reduce_dimension(sce)

print("This worked motherfuckers!!")