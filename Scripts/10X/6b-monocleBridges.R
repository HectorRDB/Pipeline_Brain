.libPaths("/home/hectorrb/.conda/envs/monocle_env/lib/R/library")
loc <- "/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_norm.rds"
output <- "/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_"

suppressMessages(library(monocle))

# Load data and convert to Çell Dataset ----
sce <- readRDS(file = loc)
pd <- new("AnnotatedDataFrame", data = as.data.frame(sce@colData))
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = rownames(sce@assays$data$counts)))
rownames(fd) <- rownames(sce@assays$data$counts)
zinbW <- sce@reducedDims[[3]]
sce <- newCellDataSet(sce@assays$data$counts,
                      phenoData = pd,
                      featureData = fd)

# Pre-process
sce <- estimateSizeFactors(sce)
sce <- estimateDispersions(sce)
print(slotNames(sce))
sce@assayData$exprs <- sce@auxOrderingData$normalize_expr_data <- t(zinbW[,1:2])
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = colnames(zinbW)[1:2]))
rownames(fd) <- colnames(zinbW)[1:2]
sce@featureData <- fd

saveRDS(sce, file = paste0(output, "monocle.rds"))
saveRDS(zinbW, file = paste0(output, "zinb.rds"))