# Packages and helper scripts ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(here)
})

source(here("Scripts", "Replicability", "meta_components.R"))
source(here("Scripts", "Replicability", "graph_visualization.R"))
source(here("Scripts", "Replicability", "01-data.R"))

# Helper functions ---- 
compute_replicability <- function(dataset, label_matrix, output_dir) {
  label_matrix <- dplyr::select(label_matrix, -dataset, -cells)
  label_sets <- colnames(label_matrix)
  
  is_not_noise <- dataset$class_label != "Noise"
  is_hvg <- rowData(dataset)$is_hvg
  
  my_dataset <- dataset[is_hvg, is_not_noise]
  my_labels <- label_matrix[is_not_noise, ]
  for (current_set in label_sets) {
    labels <- my_labels[, current_set]
    is_not_na <- !is.na(labels)
    if (sum(is_not_na) == 0) next
    stats <- analyze_components(my_dataset[, is_not_na], labels[is_not_na])
    export_components(stats, file.path(output_dir, current_set))
  }
}

analyze_components <- function(dataset, labels) {
  is_nonzero_cell <- Matrix::colSums(assay(dataset)) != 0
  dataset <- dataset[, is_nonzero_cell]
  labels <- labels[is_nonzero_cell]
  
  best_hits <- compute_best_hits(dataset, labels)
  components <- extract_components(best_hits, 0.6)
  
  return(list(
    best_hits = best_hits,
    components = components,
    labels = paste(dataset$study_id, labels, sep = "|")
  ))
}

export_components <- function(component_obj, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(component_obj$best_hits, file.path(output_dir, "best_hits.txt"))
  plot_components(component_obj$best_hits, component_obj$components$modules,
                  output_dir)
  write_component_summary(component_obj$components, output_dir)
  export_filtered_best_hits(component_obj$best_hits,
                            component_obj$components$modules, output_dir)
  
  graph <- make_directed_graph(component_obj$best_hits, 0.6, 1)
  graph <- color_graph(graph, component_obj$labels)
  pdf(file.path(output_dir, "graph_visualization.pdf"))
  plot_directed_graph(graph, 1)
  dev.off()
}

analyze_smart_tenx <- function(dataset, label_matrix, output_dir) {
  compute_replicability(dataset, label_matrix,
                        file.path(output_dir, "smart_tenx"))
}

analyze_smart <- function(dataset, label_matrix, output_dir) {
  keep <- dataset$study_id %in% c("zeng_smart_cells", "zeng_smart_nuclei")
  compute_replicability(dataset[, keep], label_matrix[keep, ],
                        file.path(output_dir, "smart"))
}

analyze_tenx <- function(dataset, label_matrix, output_dir) {
  keep <- dataset$study_id %in% c("zeng_10x_cells", "zeng_10x_nuclei")
  compute_replicability(dataset[, keep], label_matrix[keep, ],
                        file.path(output_dir, "tenx"))
}

analyze_cells <- function(dataset, label_matrix, output_dir) {
  keep <- dataset$study_id %in% c("zeng_10x_cells", "zeng_smart_cells")
  compute_replicability(dataset[, keep], label_matrix[keep, ],
                        file.path(output_dir, "cells"))
}

analyze_nuclei <- function(dataset, label_matrix, output_dir) {
  keep <- dataset$study_id %in% c("zeng_10x_nuclei", "zeng_smart_nuclei")
  compute_replicability(dataset[, keep], label_matrix[keep, ],
                        file.path(output_dir, "nuclei"))
}

# Main functions ----
## Dune ----
analyze_full_data <- function(data_path = here("data"),
                              output_dir = here("data", "Replicability",
                                                "mn_results")) {
  dataset <- load_data()
  labels <- load_labels(colnames(dataset), data_path)

  analyze_smart_tenx(dataset, labels, output_dir)
  analyze_smart(dataset, labels, output_dir)
  analyze_tenx(dataset, labels, output_dir)
  analyze_cells(dataset, labels, output_dir)
  analyze_nuclei(dataset, labels, output_dir)
}

## Hierarchical ----
analyze_all_single_merge <- function(data_path = here("data"),
                                 output_dir = here("data", "Replicability",
                                                   "mn_results")) {
  dataset <- load_smart_data()
  # Normal single Merge
  labels <- load_single_merge_labels(colnames(dataset), data_path)
  analyze_smart(dataset, labels, output_dir)
  # Comparison 2
  labels <- load_single_merge_labels(colnames(dataset), data_path,
                                     size = "large2")
  analyze_smart(dataset, labels, output_dir)
  # Comparison 3
  labels <- load_single_merge_labels(colnames(dataset), data_path,
                                     size = "large3")
  analyze_smart(dataset, labels, output_dir)
}

## single Method ----
analyze_single_methods <- function(data_path = here("data"),
                                   output_dir = here("data", "Replicability",
                                                     "mn_results")) {
  dataset <- load_smart_data()
  # Seurat
  labels <- load_single_seurat_labels(colnames(dataset), data_path)
  analyze_smart(dataset, labels, output_dir)
  # SC3
  labels <- load_single_sc3_labels(colnames(dataset), data_path)
  analyze_smart(dataset, labels, output_dir)
  # Monocle
  labels <- load_single_monocle_labels(colnames(dataset), data_path)
  analyze_smart(dataset, labels, output_dir)
}

analyze_single_methods_10x <- function(data_path = here("data"),
                                       output_dir = here("data", "Replicability",
                                                         "mn_results")) {
  dataset <- load_data()
  # Seurat
  labels <- load_seurat_all_labels(colnames(dataset), data_path)
  analyze_smart_tenx(dataset, labels, output_dir)
  # Monocle
  labels <- load_monocle_all_labels(colnames(dataset), data_path)
  analyze_smart_tenx(dataset, labels, output_dir)
}

## To run ----
main <- function() {
  analyze_single_methods()
  analyze_single_methods_10x()
  analyze_all_single_merge()
  analyze_full_data()
}

if (!interactive()) {
  main()
}