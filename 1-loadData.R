library(stringr, lib.loc = "/system/linux/lib/R-18.04/3.5/x86_64/site-library")
loc <- "/scratch/users/singlecell/MiniAtlas/data/"
# dataset <- "SMARTer_cells_MOp/"
# Can be one of SMARTer_cells_MOp/, SMARTer_nuclei_MOp/,
# 10x_cells_MOp/ and 10x_nuclei_MOp/
# Usually decided by a previous script
counts <- read.csv(paste0(loc, dataset, "exon.counts.csv.gz"),
                    header = T, row.names = 1)
colnames(counts) <- str_replace_all(colnames(counts), "\\.", "-")
meta <- read.csv(paste0(loc, dataset, "sample_metadata.csv.gz"),
                 header = T, row.names = 1)
allenClusters <- read.csv(paste0(loc, dataset, "cluster.membership.csv"),
                          header = T, col.names = c("sample", "clusters"))
meta$allenClusters <- allenClusters$clusters

library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts,
                                          logcounts = log1p(counts)),
                            colData = meta)
# saveRDS(sce, file = paste0(loc, str_replace(dataset, "/", ""), ".rds"))
# in case fails in next steps...