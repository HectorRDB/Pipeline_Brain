library(tidyverse)
loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/SMARTer_nuclei_MOp_zinbWs.rds"

# SMARTer-Nuclei ----
# Pre
# sc3    Rsec   allen seurat1 seurat2 
# 30      24      49      15      17 
# Post
clusMat <- readRDS("../../data/Smart-Seq/SMARTer_nuclei_MOp_Consensus_Clustering.rds")
out <- apply(clusMat$currentMat, 2, function(x) length(unique(x)))
df <- data.frame(methods = names(out),
                 before = c(30, 24, 49, 15, 17),
                 after = out) %>%
      gather(key = "time", value = "Nb", -methods) %>%
      mutate(time = factor(time, levels = c("before", "after")))
ggplot(df, aes(x = methods, y = Nb, fill = time)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_viridis_d(option = "E") +
  labs(x = "Clustering Methods", y = "Number of clusters", fill = "") +
  ggtitle("Reduction in number of clusters with consensus merging")
ggsave("../../Figures/Smart-Seq/SMARTer_nuclei_MOp_clusters_nb.pdf")


# SMARTer-cell ----
# Pre
# sc3    Rsec   allen seurat1 seurat2 
# 27      40      62      13      19 
# Post
clusMat <- readRDS("../../data/Smart-Seq/SMARTer_cells_MOp_Consensus_Clustering.rds")
out <- apply(clusMat$currentMat, 2, function(x) length(unique(x)))
df <- data.frame(methods = names(out),
                 before = c(27, 40, 62, 13, 19),
                 after = out) %>%
  gather(key = "time", value = "Nb", -methods) %>%
  mutate(time = factor(time, levels = c("before", "after")))
ggplot(df, aes(x = methods, y = Nb, fill = time)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_viridis_d(option = "E") +
  labs(x = "Clustering Methods", y = "Number of clusters", fill = "") +
  ggtitle("Reduction in number of clusters with consensus merging")
ggsave("../../Figures/Smart-Seq/SMARTer_cells_MOp_clusters_nb.pdf")
