library(SingleCellExperiment)
library(dplyr)
df <- readRDS("/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_sc3.rds")
df <- data.frame(cells = colnames(df), 
                 "100" = df$sc3_100_clusters)
print(dim(df))
df2 <- read.csv("/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_nuclei_MOp_sc3.csv",
                row.names = 1)
df <- inner_join(df, df2)
print(dim(df))
print(colnames(df))
write.csv(df, "/home/hectorrb/Pipeline_Brain/data/singleMethod/10x_nuclei_MOp_sc3.csv")