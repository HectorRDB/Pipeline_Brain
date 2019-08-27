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
ggplot(rep, aes(x = (replicable_clusters + non_replicable_clusters) / 2,
                y = fraction_replicable_cells, col = type)) +
  geom_point(size = 2) +
  my_theme() +
  labs(x = "Number of clusters", title = "Comparing reproducibility rates")
ggsave(here::here("Figures", "Replicability", "overall.pdf"))

ggplot(rep %>% filter(level == "Initial"),
       aes(x = (replicable_clusters + non_replicable_clusters) / 2,
                y = fraction_replicable_cells, col = type,
                shape = clustering_method)) +
  geom_point(size = 2) +
  my_theme() +
  labs(x = "Number of clusters", title = "Comparing reproducibility rates")
ggsave(here::here("Figures", "Replicability", "initial.pdf"))
