library(SingleCellExperiment)
df <- readRDS("/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_sc3.rds")
print(colData(df))
print(colnames(df)[1:10])