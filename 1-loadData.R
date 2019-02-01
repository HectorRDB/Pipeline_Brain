library(stringr)
loc <- "/scratch/users/singlecell/MiniAtlas/data/"
dataset <- "SMARTer_cells_MOp/"
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

allenMeta <- read.csv(paste0(loc, dataset, "cluster.annotation.csv"),
                      header = T, row.names = 1)

library(dplyr)

allenMetaInhibit <- allenMeta %>% filter(class_label == "GABAergic") 
inhibit <- allenClusters$clusters %in% allenMetaInhibit$cluster_id
allenMetaExcite <- allenMeta %>% filter(class_label == "Glutamatergic")
excite <- allenClusters$clusters %in% allenMetaExcite$cluster_id

library(SingleCellExperiment)
sceInhibit <- SingleCellExperiment(assays = list(counts = counts[, inhibit] %>%
                                                 as.matrix(),
                                          logcounts = as.matrix(log1p(
                                                           counts[, inhibit]))),
                                   colData = meta[inhibit, ])

sceEexcite <- SingleCellExperiment(assays = list(counts = counts[, excite] %>%
                                                   as.matrix(),
                                                 logcounts = as.matrix(log1p(
                                                   counts[, excite]))),
                                   colData = meta[excite, ])
sce <- list(Inhibit = sceInhibit,
            Excite = sceEexcite)

saveRDS(sce, file = paste0(loc, "rds/", str_replace(dataset, "/", ""), ".rds"))
# in case fails in next steps...

rm(list = setdiff(ls(), c("sce", "loc", "dataset")))