suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

data <- readRDS(here("data", "full_data.rds"))
data <- data[, data$study_id != "macosko_10x_nuclei_v3"]
df <- colData(data) %>%
  as.data.frame() %>%
  mutate(cells = colnames(data),
    study_id = case_when(study_id == "zeng_10x_cells_v3" ~ "zeng_10x_v3_cells",
                         study_id == "zeng_10x_nuclei_v3" ~ "zeng_10x_v3_nuclei",
                         study_id == "zeng_10x_cells_v2" ~ "zeng_10x_cells",
                         study_id == "zeng_10x_nuclei_v2" ~ "zeng_10x_nuclei",
                         TRUE ~ study_id),
    samples = case_when(str_detect(study_id, "v3") ~ str_extract(cells,
                                                                 "^[A-Z]*-[0-9]*"),
                        TRUE ~ cells))
rownames(df) <- df$cells
df <- df[colnames(data), ]
colnames(data) <- df$samples
data$study_id <- df$study_id
rowData(data) <- DataFrame(is_hvg = rep(TRUE, nrow(data)))
saveRDS(data, file = here("data", "full_data.rds"))