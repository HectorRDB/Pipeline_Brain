library(tidyverse)
library(here)
library(DailyHRB)

folders <- c("cells", "nuclei", "tenx", "smart")
rep <- lapply(folders, function(type) {
  df <- read.table(here::here("Replicability", type,
                              "consensus_cluster_replicability.txt"))
  df <- df %>% 
    # filter(level == "Initial") %>%
    mutate(type = type,
           clustering_method = str_remove(clustering_method, "\\..+"))
  return(df)
})
rep <- do.call("rbind", rep)
rep <- rep %>%
  mutate(Comparison = case_when(type == "cells" ~ "Between cell datasets",
                                type == "nuclei" ~ "Between nuclei datasets",
                                type == "smart" ~ "Between Smart-Seq datasets",
                                type == "tenx" ~ "Between 10x datasets"))
ggplot(rep, aes(x = (replicable_clusters + non_replicable_clusters) / 2,
                y = fraction_replicable_cells)) +
  geom_point(size = 4, alpha = .8, aes(col = Comparison)) +
  my_theme() +
  labs(x = "Number of clusters", title = "Comparing reproducibility rates",
       y = "Replicability")
ggsave(here::here("Figures", "Replicability", "overall.pdf"))

ggplot(rep %>% filter(level == "Initial"),
       aes(x = (replicable_clusters + non_replicable_clusters) / 2,
                y = fraction_replicable_cells, col = Comparison,
                shape = clustering_method)) +
  geom_point(size = 4) +
  my_theme() +
  labs(x = "Number of clusters", title = "Comparing reproducibility rates",
       y = "Replicability")
ggsave(here::here("Figures", "Replicability", "initial.pdf"))
