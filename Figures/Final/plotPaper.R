library(tidyverse)
library(here)
library(ggrepel)
isPareto <- function(x, df) {
  return(
    !any(df$fraction_replicable_cells > x[7] &
        df$nb_clusters > x[10]
    )
  )
}
df <- read.table(here("Figures", "Final", "consensus_cluster_replicability.txt"),
                 stringsAsFactors = FALSE)
df <- df %>% mutate(nb_clusters = replicable_clusters + non_replicable_clusters)
df$isPareto <- apply(df, 1, isPareto, df = df)
df$isPareto <- as.character(df$isPareto)
df$level[df$level == "Initial"] <- 0
df$level[df$level == "Final"] <- 100
df$level <- as.numeric(df$level)
df <- df %>% arrange(level)
ggplot(df, aes(x = nb_clusters,
               y = fraction_replicable_cells)) +
  geom_point(aes(color = isPareto,
                 alpha = isPareto),
             size = 4) +
  theme_classic() +
  geom_path(aes(group = clustering_name),
            alpha = .1, arrow = arrow(angle = 10,
                                      length = unit(0.15, "in"))) +
  labs(x = "Number of clusters",
       y = "fraction of cells in replicable clusters",
       col = "Pareto equilibrium\npoint") +
  scale_color_manual(values = c("grey", "red")) +
  scale_alpha_manual(values = c(.5, 1)) +
  geom_text_repel(data = df %>% filter(level == 0),
                  aes(label = clustering_name),
                  alpha = .5, size = 3) +
  guides(alpha = FALSE)
ggsave(filename = here("Figures", "Final", "figure_proposition.pdf"))
