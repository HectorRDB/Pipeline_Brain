library(here)
library(tidyverse)
# Dataset for Vignette

df <- read.csv(here("data", "Dune", "SMARTer_nuclei_MOp.csv"))
df <- df %>% select(sc3.Initial, Monocle.Initial, Seurat.Initial, cells)
allen <- read.csv(
  here("data", "Smart-Seq", "SMARTer_nuclei_MOp_cluster.membership.csv"),
  col.names = c("cells", "cluster_id"))
types <- read.csv(
  here("data", "Smart-Seq", "SMARTer_nuclei_MOp_cluster.annotation.csv"))
tsne <- read.csv(here("data", "Smart-Seq", "SMARTer_nuclei_MOp_tnse.csv")) %>%
  select(cells, x, y)
allen <- full_join(allen, types) %>% select(cells, class_label,
                                            subclass_label)
df <- inner_join(allen, df) %>%
  full_join(tsne) %>%
  filter(class_label == "GABAergic") %>%
  select(-class_label) %>%
  rename("SC3" = sc3.Initial,
         "Monocle" = Monocle.Initial,
         "Seurat" = Seurat.Initial)
save(df, file = "../Dune/data/nuclei.rda")
