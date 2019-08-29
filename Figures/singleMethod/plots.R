library(tidyverse)
library(here)
library(DailyHRB)

# Load and clean the results from doing single mergings ----
singleMethod <- read.table(here("Replicability", "singleMethod", "smart",
                               "consensus_cluster_replicability.txt")) %>%
  mutate(
    # level = str_extract(clustering_method, "\\..+$"),
    # level = str_remove(level, "^\\."),
    # level = if_else(nchar(level) == 1, paste0(level, "0"), level),
    # level = as.numeric(level),
    # clustering_method = str_remove(clustering_method, "\\..+$"),
    # clustering_method = case_when(clustering_method == "seurat" ~ "Seurat",
    #                               clustering_method == "Rsec" ~ "RSEC",
    #                               TRUE ~ clustering_method),
    merging_type = "singleMethod")

# Load and clean the results from ARI mergings
toRank <- function(i) {
  case_when(
    i == "Initial" ~ 1,
    i == 33 ~ 2,
    i == 66 ~ 3,
    i == 90 ~ 4,
    i == "Final" ~ 5
  )
}

ARIMerge <- read.table(here("Replicability", "smart",
                               "consensus_cluster_replicability.txt")) %>%
  mutate(
    # level = str_extract(clustering_method, "\\..+$"),
    # level = str_remove(level, "^\\."),
    # level = toRank(level),
    clustering_method = str_remove(clustering_method, "\\..+$"),
    merging_type = "ARI") %>%
  filter(clustering_method != "Consensus")

Init <- ARIMerge %>%
  filter(level == 1) %>%
  mutate(merging_type = "hierarchical",
         level = 0)

# singleMethod ----
## Seurat
### K.Param
ggplot(singleMethod %>% 
         mutate(level = 
                  str_extract(clustering_method, ",.*$") %>%
                  str_remove(., "^,") %>%
                  factor(., levels = c("30", "50", "100")),
                clustering_method = 
                  str_extract(clustering_method, "^.*,") %>%
                  str_remove(., "Seurat\\.") %>%
                  str_remove(., ",$") %>%
                  str_replace(., "_", ".")) %>%
         arrange(level),
       aes(
         x = (replicable_clusters + non_replicable_clusters) / 2,
         y = fraction_replicable_cells,
         col = clustering_method, shape = level)) +
  geom_point() +
  geom_path(aes(group = factor(clustering_method))) +
  my_theme() +
  labs(x = "# of clusters", y = "Replicability", col = "resolution",
       shape = "k param", title = "Seurat, varying the k parameter")
ggsave(here("Figures", "singleMethod", "singleSeurat_KParam.pdf"))

### RESOLUTION
ggplot(singleMethod %>% 
         mutate(level = 
                  str_extract(clustering_method, ",.*$") %>%
                  str_remove(., "^,") %>%
                  factor(., levels = c("30", "50", "100")),
                clustering_method = 
                  str_extract(clustering_method, "^.*,") %>%
                  str_remove(., "Seurat\\.") %>%
                  str_remove(., ",$") %>%
                  str_replace(., "_", ".")) %>%
         arrange(clustering_method),
       aes(
         x = (replicable_clusters + non_replicable_clusters) / 2,
         y = fraction_replicable_cells,
         col = level)) +
  geom_point() +
  geom_path(aes(group = factor(level))) +
  my_theme() +
  labs(x = "# of clusters", y = "Replicability", col = "k param",
       title = "Seurat, varying the resolution parameter")
ggsave(here("Figures", "singleMethod", "singleSeurat_resolution.pdf"))

# singleMethod versus ARI-merging ----
## Seurat 
# Plotting
df <- bind_rows(singleMethod, ARIMerge %>% filter(clustering_method == "Seurat")#, Init)
) %>%
  arrange(merging_type, clustering_method#, level)
  )

ggplot(df, aes(x = (replicable_clusters + non_replicable_clusters) / 2,
               y = fraction_replicable_cells,
               col = merging_type)) +
  geom_point() +
  my_theme() +
  scale_color_brewer(type = "qual") +
  labs(x = "# of clusters", y = "Replicability", col = "Type of\nmerging") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))