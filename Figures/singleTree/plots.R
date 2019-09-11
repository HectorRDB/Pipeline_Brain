library(tidyverse)
library(here)
library(DailyHRB)

# Load and clean the results from doing single mergings 
singleMerge <- read.table(here("Replicability", "singleMerge", "smart",
                               "consensus_cluster_replicability.txt")) %>%
  mutate(level = str_extract(clustering_method, "\\..+$"),
         level = str_remove(level, "^\\."),
         level = if_else(nchar(level) == 1, paste0(level, "0"), level),
         level = as.numeric(level),
         clustering_method = str_remove(clustering_method, "\\..+$"),
         clustering_method = case_when(clustering_method == "seurat" ~ "Seurat",
                                       clustering_method == "Rsec" ~ "RSEC",
                                       TRUE ~ clustering_method),
         merging_type = "hierarchical")

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
  mutate(level = str_extract(clustering_method, "\\..+$"),
         level = str_remove(level, "^\\."),
         level = toRank(level),
         clustering_method = str_remove(clustering_method, "\\..+$"),
         merging_type = "ARI") %>%
  filter(clustering_method != "Consensus")

Init <- ARIMerge %>%
  filter(level == 1, clustering_method != "RSEC") %>%
  mutate(merging_type = "hierarchical",
         level = 0)
  


# Plotting
df <- bind_rows(singleMerge, ARIMerge, Init) %>%
  arrange(merging_type, clustering_method, level)

ggplot(df, aes(x = (replicable_clusters + non_replicable_clusters) / 2,
                y = fraction_replicable_cells, col = clustering_method,
               linetype = merging_type,
               group = interaction(clustering_method, merging_type),
               label = level)) +
  geom_path(size = 3) +
  # geom_text() +
  my_theme() +
  scale_color_brewer(type = "qual") +
  labs(x = "# of clusters", y = "Replicability", col = "Clustering\nmethod",
       linetype = "Type of\nmerging") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)) 
ggsave(filename = here("Figures", "singleMerge", "singleMergeVersusARIMerge.pdf"))
