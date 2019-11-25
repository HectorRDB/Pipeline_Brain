library(SingleCellExperiment)
library(dplyr)
df <- readRDS("/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_sc3.rds")
print(colnames(df))
df <- data.frame(
  cells = colnames(df),
  "X100" = df$sc3_100_clusters
) %>%
  arrange(cells)
write.csv(df, "/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_cells_MOp_sc3.csv")