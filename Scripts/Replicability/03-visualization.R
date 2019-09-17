# Packages and helper scripts ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(here)
})

source(here("Scripts", "Replicability", "meta_components.R"))
source(here("Scripts", "Replicability", "01-data.R"))

# Helper functions ----
create_summary_figures <- function(label_matrix, result_path, output_dir,
                                   n_datasets) {
  # Create output dir if necessary
  if (!file.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  # add dataset prefix to labels
  label_matrix <- label_matrix %>%
    mutate_at(vars(-dataset, -cells), function(c) paste(label_matrix$dataset, c, sep = "|")) %>%
    as.data.frame()

  stats <- get_results(result_path, label_matrix)
  results <- stats$results
  score <- stats$score

  # if (exclude_rsec) {
  #   score <- score %>% filter(!startsWith(as.character(clustering_method), "RSEC"))
  #   results <- results %>% filter(!startsWith(as.character(clustering_method), "RSEC"))
  # }

  write.table(results, file.path(output_dir, "consensus_cluster_replicability.txt"))

  ggsave(file.path(output_dir, "fraction_mapping_clusters.pdf"),
         plot_fraction_mapped(results))
  ggsave(file.path(output_dir, "number_mapping_clusters.pdf"),
         plot_number_mapped(results, n_datasets))
  ggsave(file.path(output_dir, "fraction_mapping_cells.pdf"),
         plot_mapped_cells(results))
  ggsave(file.path(output_dir, "mapping_quality.pdf"),
         plot_mapping_quality(score))
  ggsave(file.path(output_dir, "fraction_replicable_vs_number_clusters.pdf"),
         plot_clusters_vs_mapped(results))

  pdf(file.path(output_dir, "merging_process.pdf"))
  grid::grid.newpage()
  grid::grid.draw(
    rbind(ggplotGrob(plot_number_total_2(results, n_datasets)),
          ggplotGrob(plot_mapped_cells_2(results)), size = "last")
  )
  dev.off()
}

get_results <- function(result_path, label_matrix) {
  method_labels <- get_method_labels(result_path)
  cluster_size <- lapply(dplyr::select(label_matrix, -cells), table)

  mapped <- c()
  mapped_cells <- c()
  unmapped_cells <- c()
  outliers <- c()
  mapping_score <- c()
  score <- list()
  for (l in method_labels) {
    components <- read_components(file.path(result_path, l))
    best_hits <- read_best_hits(file.path(result_path, l))
    mapped <- c(mapped, number_mapped_labels(components))
    mapped_cells <- c(mapped_cells,
                      number_mapped_cells(components, cluster_size[[l]]))
    unmapped_cells <- c(unmapped_cells,
                        number_unmapped_cells(components, cluster_size[[l]]))
    outliers <- c(outliers, number_outliers(components))
    score[[l]] <- score_components(best_hits,
                                   components$components, full_name = TRUE)$score
  }

  score <- data.frame(
    score = unlist(score),
    clustering_method = rep(names(score), sapply(score, length))
  )

  score <- score %>%
    mutate(clustering_name = sapply(strsplit(as.character(clustering_method),
                                             split = ".", fixed = TRUE), "[", 1)) %>%
    mutate(level = sapply(strsplit(as.character(clustering_method),
                                   split = ".", fixed = TRUE), "[", 2))
  if (any(score$level %in% c("Initial", "Final"))) {
    score$level <- factor(score$level, levels = c("Initial", "33", "66", "90", "Final"))
  } else {
    score$level <- factor(score$level)
  }

  results <- data.frame(
    clustering_method = method_labels,
    replicable_clusters = mapped,
    non_replicable_clusters = outliers,
    replicable_cells = mapped_cells,
    non_replicable_cells = unmapped_cells
  )

  results <- results %>%
    mutate(fraction_replicable_clusters = replicable_clusters / 
             (replicable_clusters + non_replicable_clusters)) %>%
    mutate(fraction_replicable_cells = replicable_cells / 
             (replicable_cells + non_replicable_cells)) %>%
    mutate(clustering_name = sapply(strsplit(as.character(clustering_method),
                                             split = ".", fixed = TRUE), "[", 1)) %>%
    mutate(level = sapply(strsplit(as.character(clustering_method),
                                   split = ".", fixed = TRUE), "[", 2))

  if (any(results$level %in% c("Initial", "Final"))) {
    results$level <- factor(results$level,
                            levels = c("Initial", "33", "66", "90", "Final"))
  } else {
    results$level <- factor(results$level)
  }

  return(list(score = score, results = results))
}

get_method_labels <- function(result_path) {
  labels <- list.dirs(result_path, full.names = FALSE)
  labels <- labels[labels != ""]
  return(labels)
}

read_components <- function(result_path) {
  components <- read.table(file.path(result_path, "component_summary.txt"),
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  comp <- lapply(unique(components$Component), function(c) {
    components$Label[components$Component == c]
  })
  names(comp) <- unique(components$Component)
  return(list(
    components = comp[names(comp) != "-1"],
    outliers = comp[["-1"]]
  ))
}

number_mapped_labels <- function(components) {
  return(length(unlist(components$components)))
}

number_mapped_cells <- function(components, cluster_sizes) {
  return(sum(cluster_sizes[unlist(components$components)]))
}

number_outliers <- function(components) {
  if (is.null(components$outliers)) {
    return(0)
  } else {
    return(sum(get_cell_type(components$outliers) != "-1"))  
  }
}

number_unmapped_cells <- function(components, cluster_sizes) {
  unmapped_labels <- components$outliers
  if (is.null(unmapped_labels)) {
    return(0)
  }
  unmapped_labels <- unmapped_labels[get_cell_type(unmapped_labels) != "-1"]
  return(sum(cluster_sizes[unmapped_labels]))
}

read_best_hits <- function(result_path) {
  result <- read.table(file.path(result_path, "best_hits.txt"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  return(as.matrix(result))
}

plot_all <- function(results, score, n_datasets) {
  print(plot_fraction_mapped(results))
  print(plot_number_mapped(results, n_datasets))
  print(plot_number_total(results, n_datasets))
  print(plot_mapped_cells(results))
  print(plot_mapping_quality(score))
  print(plot_clusters_vs_mapped(results))
}

plot_fraction_mapped <- function(results) {
  ggplot(results, aes(x = level, y = fraction_replicable_clusters,
                      group = clustering_name, col = clustering_name)) +
    geom_line() +
    geom_point() +
    theme_classic() +
    theme(legend.position = "top", axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab("Merging step") +
    ylab("Fraction of clusters replicable in at least one other dataset") +
    labs(col = "Clustering method")
}

plot_number_mapped <- function(results, n_datasets) {
  ggplot(results, aes(x = level, y = replicable_clusters / n_datasets,
                      group = clustering_name, col = clustering_name)) +
    geom_line() +
    geom_point() +
    theme_classic() +
    theme(legend.position = "top", axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab("Merging step") +
    ylab("Average number of clusters replicable in at least one other dataset") +
    labs(col = "Clustering method")
}

plot_mapped_cells <- function(results) {
  ggplot(results, aes(
    x = level, y = fraction_replicable_cells,
    group = clustering_name, col = clustering_name
  )) +
    geom_line() +
    geom_point() +
    theme_classic() +
    theme(legend.position = "top", axis.ticks.x = element_blank(), 
          axis.text.x = element_blank()) +
    xlab("Merging step") +
    ylab("Fraction of cells in replicable clusters") +
    labs(col = "Clustering method")
}

plot_number_total <- function(results, n_datasets) {
  ggplot(results, aes(
    x = level, y = (replicable_clusters + non_replicable_clusters) / n_datasets,
    group = clustering_name, col = clustering_name
  )) +
    geom_line() +
    geom_point() +
    theme_classic() +
    theme(legend.position = "top", axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab("Merging step") +
    ylab("Average number of clusters") +
    labs(col = "Clustering method")
}

plot_mapped_cells_2 <- function(results) {
  ggplot(results, aes(
    x = level, y = fraction_replicable_cells,
    group = clustering_name, label = as.numeric(level), col = clustering_name
  )) +
    geom_line(alpha = 0.5) +
    geom_text(size = 6, col = "grey30") +
    geom_text(
      data = subset(results, as.numeric(level) == 1),
      aes(label = clustering_name),
      nudge_x = -0.31
    ) +
    theme_classic() +
    theme(legend.position = "top", axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab("Merging step") +
    ylab("Fraction of cells in replicable clusters") +
    guides(col = FALSE)
}

plot_number_total_2 <- function(results, n_datasets) {
  ggplot(results, aes(
    x = level, y = (replicable_clusters + non_replicable_clusters) / n_datasets,
    group = clustering_name, label = as.numeric(level), col = clustering_name
  )) +
    geom_line(alpha = 0.5) +
    geom_text(size = 6, col = "grey30") +
    geom_text(
      data = subset(results, as.numeric(level) == 1),
      aes(label = clustering_name),
      nudge_x = -0.31
    ) +
    theme_classic() +
    theme(
      legend.position = "top", axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ylab("Average number of clusters") +
    guides(col = FALSE)
}

plot_mapping_quality <- function(score) {
  ggplot(score, aes(x = level, y = score, fill = clustering_name)) +
    geom_boxplot() +
    theme_classic() +
    xlab("Consensus level") +
    ylab("Mapping quality (average AUROC)")
}

plot_clusters_vs_mapped <- function(results) {
  results %>%
    arrange(level) %>%
    ggplot(aes(
      x = replicable_clusters + non_replicable_clusters,
      y = fraction_replicable_cells,
      col = clustering_name,
      group = clustering_name,
      shape = level
    )) +
    geom_path() +
    geom_point(size = 2) +
    theme_classic() +
    # theme(legend.position = "top") +
    xlab("Total number of clusters") +
    ylab("Fraction of cells in a replicable cluster") +
    labs(col = "Clustering method", shape = "Consensus level")
}

# Main functions ----
main_full_data <- function(
  result_path = here("data", "Replicability", "mn_results", "Dune"),
  output_dir = here("data", "Replicability", "Dune")) 
  {
    
  label_matrix <- load_labels(load_qc_cells())
  
  create_summary_figures(
    label_matrix, file.path(result_path, "smart_tenx"),
    file.path(output_dir, "smart_tenx"), 4
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "tenx"),
    file.path(output_dir, "tenx"), 2
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "cells"),
    file.path(output_dir, "cells"), 2
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "nuclei"),
    file.path(output_dir, "nuclei"), 2
  )
}

main_all_Dunes <- function(
  result_path = here("data", "Replicability", "mn_results", "Dune_Smart"),
  output_dir = here("data", "Replicability", "Dune_Smart")) 
  {
  
  dataset <- load_smart_data()
  
  # Dune normal
  label_matrix <- load_Dune_labels(load_qc_cells(here("data", "qc_cells_smart.txt")))
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix,
                         file.path(result_path, "Normal", "smart"),
                         file.path(output_dir, "Normal"), 2
  )
  # Dune large 2
  label_matrix <- load_Dune_labels(load_qc_cells(here("data", "qc_cells_smart.txt")),
                                   size = "large2")
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix,
                         file.path(result_path, "large2", "smart"),
                         file.path(output_dir, "Large2"), 2
  )
  # Dune large 3
  label_matrix <- load_Dune_labels(load_qc_cells(here("data", "qc_cells_smart.txt")),
                                   size = "large3")
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix,
                         file.path(result_path, "large3", "smart"),
                         file.path(output_dir, "Large3"), 2
  )
}

main_single_merge <- function(
  result_path = here("data", "Replicability", "mn_results", "singleTree"),
  output_dir = here("data", "Replicability", "SingleTree")) 
  {
  
  dataset <- load_smart_data()
  
  # Normal hierarchical
  label_matrix <- load_single_merge_labels(load_qc_cells(here("data", "qc_cells_smart.txt")))
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix,
                         file.path(result_path, "Normal", "smart"),
                         file.path(output_dir, "Normal"), 2
  )
  # Hierarchical large 2
  label_matrix <- load_single_merge_labels(load_qc_cells(here("data", "qc_cells_smart.txt")),
                                           size = "_large2")
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix,
                         file.path(result_path, "large2", "smart"),
                         file.path(output_dir, "Large2"), 2
  )
  # Hierarchical large 3
  label_matrix <- load_single_merge_labels(load_qc_cells(here("data", "qc_cells_smart.txt")),
                                           size = "_large3")
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix,
                         file.path(result_path, "large3", "smart"),
                         file.path(output_dir, "Large3"), 2
  )
}

main_single_method_all <- function(
  result_path = here("data", "Replicability", "mn_results", "SingleMethod"),
  output_dir = here("data", "Replicability", "SingleMethod")) 
  {
  
  dataset <- load_smart_data()
  
  # Smart-Seq only (3 methods)
  label_matrix <- inner_join(
    load_single_seurat_labels(load_qc_cells(here("data", "qc_cells_smart.txt"))),
    load_single_sc3_labels(load_qc_cells(here("data", "qc_cells_smart.txt")))
  ) %>%
    inner_join(load_single_monocle_labels(
      load_qc_cells(here("data", "qc_cells_smart.txt"))))
  label_matrix <- label_matrix[dataset$class_label != "Noise", ]
  create_summary_figures(label_matrix, file.path(result_path, "smart"),
                         file.path(output_dir, "smart"), 2)
  
  # Smart-Seq and 10x (2 methods)
  label_matrix <- inner_join(
    load_seurat_all_labels(load_qc_cells()),
    load_monocle_all_labels(load_qc_cells())
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "smart_tenx"),
    file.path(output_dir, "smart_tenx"), 4
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "tenx"),
    file.path(output_dir, "tenx"), 2
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "cells"),
    file.path(output_dir, "cells"), 2
  )
  create_summary_figures(
    label_matrix, file.path(result_path, "nuclei"),
    file.path(output_dir, "nuclei"), 2
  )
}

main <- function() {
  print("full data")
  main_full_data()
  print("all Dunes")
  main_all_Dunes()
  print("single method")
  main_single_method_all()
  print("single merges")
  main_single_merge()
}

if (!interactive()) {
  main()
}