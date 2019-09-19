suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})


main <- function() {
  create_replicability_tsnes()
}

create_replicability_tsnes <- function(
  result_path = here("data", "Replicability", "mn_results", "Dune", "smart_tenx"),
  output_dir = here("Figures/tSNE"),
  tSNE_path = here("data", "tSNE"),
  data_path = here("data"))
  {
  
  if (!file.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  labels <- load_labels(data_path)
  tsne_coord <- load_tsne_coord(labels$cells, tSNE_path)
  replicability_matrix <- compute_replicability_matrix(labels, result_path)
  plot_tsne_each_step(replicability_matrix, tsne_coord, output_dir)
  plot_tsne_global(replicability_matrix, tsne_coord, output_dir)
}

load_labels <- function(data_path) {
  data_path <- file.path(data_path, "Dune")
  label_matrix <- bind_rows(
    zeng_smart_cells = read.csv(file.path(data_path, "SMARTer_cells_MOp.csv")),
    zeng_smart_nuclei = read.csv(file.path(data_path, "SMARTer_nuclei_MOp.csv")),
    zeng_10x_cells = read.csv(file.path(data_path, "10x_cells_MOp.csv")),
    zeng_10x_nuclei = read.csv(file.path(data_path, "10x_nuclei_MOp.csv")),
    .id = "dataset"
  ) %>% select(-X)

  # add dataset prefix
  label_matrix <- label_matrix %>%
    mutate_at(vars(-dataset, -cells), function(c) paste(label_matrix$dataset, c,
                                                        sep = "|"))
  label_matrix <- as.data.frame(label_matrix)

  return(label_matrix)
}

load_tsne_coord <- function(cell_names, tSNE_path) {
  tsne_coord <- bind_rows(
    zeng_smart_cells = read.csv(
      file.path(tSNE_path, "SMARTer_cells_MOp_tnse.csv")),
    zeng_smart_nuclei = read.csv(
      file.path(tSNE_path, "SMARTer_nuclei_MOp_tnse.csv")),
    zeng_10x_cells = read.csv(
      file.path(tSNE_path, "10x_cells_MOp_tnse.csv")),
    zeng_10x_nuclei = read.csv(
      file.path(tSNE_path, "10x_nuclei_MOp_tnse.csv")),
    .id = "dataset"
  )

  tsne_coord <- tsne_coord %>%
    dplyr::select(-X)

  return(tsne_coord)
}

compute_replicability_matrix <- function(label_matrix, result_path) {
  labels <- get_method_labels(result_path)
  is_replicable <- lapply(set_names(labels), function(l) {
    find_replicable_cluster(file.path(result_path, l))
  }
  )

  replicability_matrix <- label_matrix
  for (name in names(is_replicable)) {
    replicability_matrix[[name]] <- is_replicable[[name]][replicability_matrix[[name]]]
  }

  replicability_matrix <- replicability_matrix %>%
    gather("method.step", "is_replicable", -dataset, -cells) %>%
    mutate(method = my_word(method.step, 1, ".")) %>%
    mutate(merging_step = my_word(method.step, 2, "."))

  return(replicability_matrix)
}

get_method_labels <- function(output_dir) {
  labels <- list.dirs(output_dir, full.names = FALSE)
  labels <- labels[labels != ""]
  labels <- labels[!startsWith(labels, "Allen")]
  return(labels)
}

find_replicable_cluster <- function(component_dir) {
  read_component_summary(component_dir) %>% deframe() > 0
}

read_component_summary <- function(results_dir) {
  read.table(file.path(results_dir, "component_summary.txt"), sep = "\t",
             header = TRUE, stringsAsFactors = FALSE)
}

my_word <- function(x, position, split, fixed = TRUE) {
  sapply(strsplit(x, split = split, fixed = fixed), "[", position)
}

plot_tsne_each_step <- function(replicability_matrix, tsne_coord, output_dir) {
  dataset_names <- unique(replicability_matrix$dataset)
  merging_steps <- unique(replicability_matrix$merging_step)

  for (name in dataset_names) {
    dot_size <- get_dot_size(name)
    for (step in merging_steps) {
      my_plot <- replicability_matrix %>%
        compute_scores(name, step) %>%
        plot_scores(tsne_coord, name, step, dot_size)
      ggsave(file.path(output_dir, paste0(name, "_", step, ".png")), my_plot)
    }
  }
}

get_dot_size <- function(dataset_name) {
  if (dataset_name %in% c("zeng_smart_cells", "zeng_smart_nuclei")) {
    return(1)
  } else {
    return(0.01)
  }
}

compute_scores <- function(replicability_matrix, dataset_name, merging_step_name) {
  replicability_matrix %>%
    filter(method != "Consensus") %>%
    drop_na() %>%
    filter(dataset == dataset_name & merging_step == merging_step_name) %>%
    group_by(cells) %>%
    summarize(dataset = first(dataset), replicability_score = sum(is_replicable))
}

plot_scores <- function(scores, tsne_coord, dataset_name, merging_step_name,
                        dot_size = 0.01) {
  scores %>%
    inner_join(tsne_coord) %>%
    ggplot(aes(x = x, y = y, col = replicability_score)) +
    geom_point(size = dot_size, alpha = 0.5) +
    xlab("tSNE1") +
    ylab("tSNE2") +
    theme_classic() +
    ggtitle(paste0(dataset_name, " (merging step: ", merging_step_name, ")"))
}

plot_tsne_global <- function(replicability_matrix, tsne_coord, output_dir) {
  dataset_names <- unique(replicability_matrix$dataset)

  for (name in dataset_names) {
    dot_size <- get_dot_size(name)
    global_scores <- replicability_matrix %>%
      filter(method != "Consensus") %>%
      drop_na() %>%
      filter(dataset == name) %>%
      group_by(cells) %>%
      summarize(dataset = first(dataset),
                replicability_score = sum(is_replicable))
    my_plot <- plot_scores(global_scores, tsne_coord, name, "all", dot_size)
    ggsave(file.path(output_dir, paste0(name, "_all_steps.png")), my_plot)
  }
}

if (!interactive()) {
  main()
}