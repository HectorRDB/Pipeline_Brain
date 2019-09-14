library(here)
library(tidyverse)
# Tool dataset

set.seed(24)
init <- sample(1:10, 100, replace = TRUE)
clusMat <- matrix(rep(init, 5), ncol = 5, byrow = FALSE)
clusMat[clusMat[, 2] == 8, 2] <- sample(11:12, sum(clusMat[, 2] == 8), TRUE)
clusMat[clusMat[, 3] == 8, 3] <- sample(11:12, sum(clusMat[, 3] == 8), TRUE)
clusMat[clusMat[, 2] == 4, 4] <- sample(11:13, sum(clusMat[, 4] == 4), TRUE)
clusMat[clusMat[, 2] == 2, 5] <- sample(11:14, sum(clusMat[, 5] == 2), TRUE)

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
