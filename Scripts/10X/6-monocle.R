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

suppressMessages(library(monocle))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zinbwave))

# Load data and convert to Delayed Array ----
sce <- readRDS(file = loc)
pd <- new("AnnotatedDataFrame", data = as.data.frame(sce@colData))
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(assays(sce)$counts)))
zinbW <- reducedDim(sce, type = reducedDimNames(sce)[3])
rownames(fd) <- rownames(assays(sce)$counts)
sce <- newCellDataSet(assays(sce)$counts,
                      phenoData = pd,
                      featureData = fd)

# Pre-process
sce <- estimateSizeFactors(sce)
sce <- estimateDispersions(sce)
sce@normalized_data_projection <- zinbW
sce@assayData$exprs <- sce@auxOrderingData$normalize_expr_data <- t(zinbW[,1:2])
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = colnames(zinbW)[1:2]))
rownames(fd) <- colnames(zinbW)[1:2]
sce@featureData <- fd

saveRDS(sce, file = output)
