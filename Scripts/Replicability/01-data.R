suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

# source("variable_genes.R")
# source("identifier_conversion.R")
source(here("Scripts", "Replicability", "datasets.R"))

# Helper functions ----
export_qc_cells <- function(dataset = load_data(),
                            filename = here("data", "qc_cells.txt")) {
  write(colnames(dataset[, dataset$class_label != "Noise"]), filename)
}

load_data <- function() {
  readRDS(here("data", "full_data.rds"))
}

load_qc_cells <- function(filename = here("data", "qc_cells.txt")) {
  scan(filename, "character")
}

load_lab_data <- function() {
  dataset <- readRDS(here("data", "lab_data.rds"))
  return(dataset)
}

# Load dune data ----
load_smart_data <- function() {
  result <- readRDS(here("data", "full_data.rds"))
  result <- result[,
                result$study_id %in% c("zeng_smart_cells", "zeng_smart_nuclei")]
}

load_labels <- function(cell_names, data_path = here("data")) {
  input_dir <- file.path(data_path, "Dune")
  label_matrix <- bind_rows(
    zeng_smart_cells = read.csv(file.path(input_dir, "SMARTer_cells_MOp.csv")),
    zeng_smart_nuclei = read.csv(file.path(input_dir, "SMARTer_nuclei_MOp.csv")),
    zeng_10x_cells = read.csv(file.path(input_dir, "10x_cells_MOp.csv")),
    zeng_10x_nuclei = read.csv(file.path(input_dir, "10x_nuclei_MOp.csv")),
    zeng_10x_v3_cells = read.csv(file.path(input_dir, "10x_v3_cells_MOp.csv")),
    zeng_10x_v3_nuclei = read.csv(file.path(input_dir, "10x_v3_nuclei_MOp.csv")),
    .id = "dataset"
  )

  # reorder cells to match data
  row_match <- match(cell_names, label_matrix$cells)
  label_matrix <- label_matrix[row_match, ]

  return(label_matrix)
}

load_Dune_labels <- function(cell_names, data_path = here("data"),
                             size = "normal", metric = "") {
  if (size == "normal") {
    input_dir <- file.path(data_path, "Dune")
    label_matrix <- bind_rows(
      zeng_smart_cells = 
        read.csv(file.path(input_dir, paste0("SMARTer_cells_MOp", metric, ".csv"))),
      zeng_smart_nuclei = 
        read.csv(file.path(input_dir, paste0("SMARTer_nuclei_MOp", metric, ".csv"))),
      .id = "dataset"
    )
  } else {
    input_dir <- file.path(data_path, "singleTree")
    label_matrix <- bind_rows(
      zeng_smart_cells = read.csv(
        file.path(input_dir, paste0("SMARTer_cells_MOp_", size, metric, "_Dune.csv"))),
      zeng_smart_nuclei = read.csv(
        file.path(input_dir, paste0("SMARTer_nuclei_MOp_", size, metric, "_Dune.csv"))),
      .id = "dataset"
    )
  }
  
  if ("X" %in% colnames(label_matrix)) {
    label_matrix <- label_matrix %>% select(-X)
  }
  # reorder cells to match data
  row_match <- match(cell_names, label_matrix$cells)
  label_matrix <- label_matrix[row_match, ]
  
  return(label_matrix)
}

# Load hierarchical ----

load_single_merge_labels <- function(cell_names, data_path = here("data"),
                                     size = "", type = "DE") {
  input_dir <- file.path(data_path, "singleTree")
  result <- bind_rows(
    zeng_smart_cells = read.csv(
      file.path(input_dir, 
                paste0("SMARTer_cells_MOp", size , "_hierarchical_", type,
                       ".csv"))),
    zeng_smart_nuclei = read.csv(
      file.path(input_dir,
                paste0("SMARTer_nuclei_MOp", size , "_hierarchical_", type,
                       ".csv"))),
    .id = "dataset"
  )

  # restrict to steps where both datasets have at least 2 clusters
  n_clusters_cells <- apply(result[result$dataset == "zeng_smart_cells", ], 2,
                            function(x) length(table(x)))
  n_clusters_nuclei <- apply(result[result$dataset == "zeng_smart_nuclei", ], 2,
                             function(x) length(table(x)))
  keep <- n_clusters_cells > 1 & n_clusters_nuclei > 1
  keep[1] <- TRUE
  result <- result[, keep]
  result <- result %>% select(dataset, cells, colnames(result))
  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.MERGING_STEP)
  methods <- colnames(result)[-(1:2)]
  if (type == "DE") {
    method_name <- stringr::word(methods, 1, sep = stringr::fixed("."))
    method_level <- stringr::word(methods, 3, sep = stringr::fixed("."))
    method_level[is.na(method_level)] <- "00"
    colnames(result)[-(1:2)] <- paste(method_name, method_level, sep = ".")  
  } else {
    methods <- stringr::word(methods, 2, sep = stringr::fixed("."))
    colnames(result)[-(1:2)] <- methods
  }
  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]
  
  return(result)
}

# Load single Method ----
read_single_method <- function(filename) {
  result <- read.csv(filename, check.names = FALSE, row.names = 1)
  if ("cells" %in% colnames(result)) {
    result <- as_tibble(result)
  } else {
    result <- as_tibble(result, rownames = "cells")  
  }
  return(result %>% select(cells, colnames(result)))
}

load_single_seurat_labels <- function(cell_names, data_path = here("data")) {
  input_dir <- file.path(data_path, "singleMethod")
  result <- bind_rows(
    zeng_smart_cells = read_single_method(
      file.path(input_dir, "SMARTer_cells_MOp_Seurat.csv")),
    zeng_smart_nuclei = read_single_method(
      file.path(input_dir, "SMARTer_nuclei_MOp_Seurat.csv")),
    .id = "dataset"
  )

  # restrict to steps where both datasets have at least 2 clusters
  n_clusters_cells <- apply(result[result$dataset == "zeng_smart_cells", ], 2,
                            function(x) length(table(x)))
  n_clusters_nuclei <- apply(result[result$dataset == "zeng_smart_nuclei", ], 2,
                             function(x) length(table(x)))
  keep <- n_clusters_cells > 1 & n_clusters_nuclei > 1
  keep[1] <- TRUE
  result <- result[, keep]

  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.PARAMETER)
  methods <- colnames(result)[-(1:2)]
  method_name <- "Seurat"
  method_level <- gsub("\\.", "\\_", methods)
  colnames(result)[-(1:2)] <- paste(method_name, method_level, sep = ".")

  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]

  return(as.data.frame(result))
}

load_single_sc3_labels <- function(cell_names, data_path = here("data")) {
  input_dir <- file.path(data_path, "singleMethod")
  # Rename cells
  zeng_smart_cells <- read_single_method(
    file.path(input_dir, "SMARTer_cells_MOp_SC3.csv")) %>%
    mutate(dataset = "zeng_smart_cells")
  zeng_smart_cells <- zeng_smart_cells %>%
    select(.data = ., "dataset", colnames(zeng_smart_cells))
  # Rename nuclei
  zeng_smart_nuclei <- read_single_method(
    file.path(input_dir, "SMARTer_nuclei_MOp_SC3.csv")) %>%
    mutate(dataset = "zeng_smart_nuclei")
  zeng_smart_nuclei <- zeng_smart_nuclei %>%
    select(.data = ., "dataset", colnames(zeng_smart_nuclei))
  
  result <- bind_rows(
    zeng_smart_cells,
    zeng_smart_nuclei
  )

  # restrict to steps where both datasets have at least 2 clusters
  n_clusters_cells <- apply(result[result$dataset == "zeng_smart_cells", ], 2,
                            function(x) length(table(x)))
  n_clusters_nuclei <- apply(result[result$dataset == "zeng_smart_nuclei", ], 2,
                             function(x) length(table(x)))
  keep <- n_clusters_cells > 1 & n_clusters_nuclei > 1
  keep[1] <- TRUE
  result <- result[, keep]

  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.PARAMETER)
  colnames(result)[-(1:2)] <- paste("SC3", colnames(result)[-(1:2)], sep = ".")

  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]

  return(as.data.frame(result))
}

load_single_monocle_labels <- function(cell_names, data_path = here("data")) {
    input_dir <- file.path(data_path, "singleMethod")
    result <- bind_rows(
        zeng_smart_cells = read_single_method(
            file.path(input_dir, "SMARTer_cells_MOp_Monocle.csv")),
        zeng_smart_nuclei = read_single_method(
            file.path(input_dir, "SMARTer_nuclei_MOp_Monocle.csv")),
        .id = "dataset"
    )
    
    # restrict to steps where both datasets have at least 2 clusters
    n_clusters_cells <- apply(result[result$dataset == "zeng_smart_cells", ], 2,
                              function(x) length(table(x)))
    n_clusters_nuclei <- apply(result[result$dataset == "zeng_smart_nuclei", ], 2,
                               function(x) length(table(x)))
    keep <- n_clusters_cells > 1 & n_clusters_nuclei > 1
    keep[1] <- TRUE
    result <- result[, keep]
    
    # rename columns for compatibility with analysis/visualization modules
    # (format METHOD_NAME.PARAMETER)
    colnames(result)[-(1:2)] <- paste("Monocle",
                                      colnames(result)[-(1:2)], sep = ".")
    
    # reorder cells to match data
    row_match <- match(cell_names, result$cells)
    result <- result[row_match, ]
    
    return(as.data.frame(result))
}

load_seurat_all_labels <- function(cell_names, data_path = here("data")) {
  input_dir <- file.path(data_path, "singleMethod")
  result <- bind_rows(
    zeng_smart_cells = read_single_method(
      file.path(input_dir, "SMARTer_cells_MOp_Monocle.csv")),
    zeng_smart_nuclei = read_single_method(
      file.path(input_dir, "SMARTer_nuclei_MOp_Monocle.csv")),
    zeng_10x_cells = read_single_method(
      file.path(input_dir, "10x_cells_MOp_Monocle.csv")),
    zeng_10x_nuclei = read_single_method(
      file.path(input_dir, "10x_nuclei_MOp_Monocle.csv")),
    zeng_10x_v3_cells = read_single_method(
      file.path(input_dir, "10x_v3_cells_MOp_Monocle.csv")),
    zeng_10x_v3_nuclei = read_single_method(
      file.path(input_dir, "10x_v3_nuclei_MOp_Monocle.csv")),
    .id = "dataset"
  )
  
  # restrict to steps where all datasets have at least 2 clusters
  n_clusters_cells_smart <- apply(result[result$dataset == "zeng_smart_cells", ],
                                  2, function(x) length(table(x)))
  n_clusters_nuclei_smart <- apply(result[result$dataset == "zeng_smart_nuclei", ],
                                   2, function(x) length(table(x)))
  n_clusters_cells_10x <- apply(result[result$dataset == "zeng_10x_cells", ],
                                2, function(x) length(table(x)))
  n_clusters_nuclei_10x <- apply(result[result$dataset == "zeng_10x_nuclei", ],
                                 2, function(x) length(table(x)))
  n_clusters_cells_10x_v3 <- apply(result[result$dataset == "zeng_10x_v3_cells", ],
                                   2, function(x) length(table(x)))
  n_clusters_nuclei_10x_v3 <- apply(result[result$dataset == "zeng_10x_v3_nuclei", ],
                                    2, function(x) length(table(x)))
  keep <- n_clusters_cells_smart > 1 & n_clusters_nuclei_smart > 1 &
    n_clusters_cells_10x > 1 & n_clusters_nuclei_10x > 1 &
    n_clusters_cells_10x_v3 > 1 & n_clusters_nuclei_10x_v3 > 1
  keep[1] <- TRUE
  result <- result[, keep]
  
  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.PARAMETER)
  methods <- colnames(result)[-(1:2)]
  method_name <- "Seurat"
  method_level <- gsub("\\.", "\\_", methods)
  colnames(result)[-(1:2)] <- paste(method_name, method_level, sep = ".")
  
  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]
  
  return(as.data.frame(result))
}

load_monocle_all_labels <- function(cell_names, data_path = here("data")) {
  input_dir <- file.path(data_path, "singleMethod")
  result <- bind_rows(
    zeng_smart_cells = read_single_method(
      file.path(input_dir, "SMARTer_cells_MOp_Monocle.csv")),
    zeng_smart_nuclei = read_single_method(
      file.path(input_dir, "SMARTer_nuclei_MOp_Monocle.csv")),
    zeng_10x_cells = read_single_method(
      file.path(input_dir, "10x_cells_MOp_Monocle.csv")),
    zeng_10x_nuclei = read_single_method(
      file.path(input_dir, "10x_nuclei_MOp_Monocle.csv")),
    zeng_10x_v3_cells = read_single_method(
      file.path(input_dir, "10x_v3_cells_MOp_Monocle.csv")),
    zeng_10x_v3_nuclei = read_single_method(
      file.path(input_dir, "10x_v3_nuclei_MOp_Monocle.csv")),
    .id = "dataset"
  )
  
  # restrict to steps where all datasets have at least 2 clusters
  n_clusters_cells_smart <- apply(result[result$dataset == "zeng_smart_cells", ],
                                  2, function(x) length(table(x)))
  n_clusters_nuclei_smart <- apply(result[result$dataset == "zeng_smart_nuclei", ],
                                   2, function(x) length(table(x)))
  n_clusters_cells_10x <- apply(result[result$dataset == "zeng_10x_cells", ],
                                2, function(x) length(table(x)))
  n_clusters_nuclei_10x <- apply(result[result$dataset == "zeng_10x_nuclei", ],
                                 2, function(x) length(table(x)))
  n_clusters_cells_10x_v3 <- apply(result[result$dataset == "zeng_10x_v3_cells", ],
                                2, function(x) length(table(x)))
  n_clusters_nuclei_10x_v3 <- apply(result[result$dataset == "zeng_10x_v3_nuclei", ],
                                 2, function(x) length(table(x)))
  keep <- n_clusters_cells_smart > 1 & n_clusters_nuclei_smart > 1 &
    n_clusters_cells_10x > 1 & n_clusters_nuclei_10x > 1 &
    n_clusters_cells_10x_v3 > 1 & n_clusters_nuclei_10x_v3 > 1
  keep[1] <- TRUE
  result <- result[, keep]
  
  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.PARAMETER)
  colnames(result)[-(1:2)] <- paste("Monocle",
                                    colnames(result)[-(1:2)], sep = ".")
  
  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]
  
  return(as.data.frame(result))
}
# Load lab comp ----
load_Lab_labels <- function(data_path = "data") {
    input_dir <- file.path(data_path, "Dune")
    label_matrix <- bind_rows(
      Zeng = read.csv(file.path(input_dir, "10x_nuclei_MOp.csv")),
      Regev = read.csv(file.path(input_dir, "Regev.csv")),
      .id = "dataset"
    )
  
  # reorder cells to match data
  
  return(label_matrix)
}