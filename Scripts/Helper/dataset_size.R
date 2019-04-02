library(clusterExperiment)
library(ggplot2)

# Load data ----
loc <- "/pylon5/ib5phhp/hectorrb/ProcessedData/10x_nuclei_MOp_filt.rds"
sce <- readRDS(file = loc)
print(sce)

loc <- "/pylon5/ib5phhp/hectorrb/ProcessedData/10x_cells_MOp_filt.rds"
sce <- readRDS(file = loc)
print(sce)