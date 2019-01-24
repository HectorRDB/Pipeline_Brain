# dataset <- "SMARTer_cells_MOp/"
source("1-loadData.R")

counts <- data.frame(assays(sce)$counts)
counts[is.na(counts)] <- 0
assays(sce)$counts <- counts

filt <- apply(assays(sce)$counts, 1, function(x) {
  sum(x >= 50) >= 50
})

sce <- sce[filt, ]

# saveRDS(sce, file = paste0(loc, str_replace(dataset, "/", ""), "_filt.rds"))
# in case fails in next steps...