library(Rtsne)
library(SingleCellExperiment)
library(tidyverse)
loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_zinbWs.rds"
sce <- readRDS(loc)

dims <- reducedDimNames(sce)
zinbWs <- reducedDims(sce)
clusters <- read.csv(
"/scratch/users/singlecell/MiniAtlas/data/SMARTer_nuclei_MOp/cluster.annotation.csv",
header = T)
cols2 <- clusters$cluster_color %>% as.character()
names(cols2) <- clusters$cluster_id %>% as.character()

# Is zinbWave the problem? ----
for (i in 1:length(zinbWs)) {
  zinbW <- zinbWs[[i]]
  print(i)
  TNSE <- Rtsne(zinbW)
  print("     tsne")
  df <- data.frame(x = TNSE$Y[, 1], y = TNSE$Y[, 2],
                   cols = as.factor(colData(sce)$allenClusters))
  p <- ggplot(df, aes(x = x, y = y, col = cols)) +
        geom_point() +
        theme_classic() +
        scale_color_manual(values = cols2, breaks = names(cols2)) +
        labs(x = "dim1", y = "dim2")
  ggsave(paste0(
  "/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/Figures/Exploration/tsne_nuclei_K", dims[i], ".pdf"),
  p)
  print("     saving plot")
}
