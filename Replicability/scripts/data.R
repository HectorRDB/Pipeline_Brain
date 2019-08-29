suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(tidyverse)
})


export_qc_cells = function(dataset = load_data(), filename = "qc_cells.txt") {
    write(colnames(dataset[, dataset$class_label != "Noise"]), filename)
}

load_data = function() {
    readRDS("full_data.rds")
}

load_qc_cells = function(filename = "qc_cells.txt") {
    scan(filename, "character")
}

load_smart_data = function() {
    result = readRDS("full_data.rds")
    result = result[, result$study_id %in% c("zeng_smart_cells", "zeng_smart_nuclei")]
}

load_labels = function(cell_names, data_path = "../../data") {
    input_dir = file.path(data_path, "ClusterLabels")
    label_matrix = bind_rows(
        zeng_smart_cells = read.csv(file.path(input_dir, "SMARTer_cells_MOp.csv")),
        zeng_smart_nuclei = read.csv(file.path(input_dir, "SMARTer_nuclei_MOp.csv")),
        zeng_10x_cells = read.csv(file.path(input_dir, "10x_cells_MOp.csv")),
        zeng_10x_nuclei = read.csv(file.path(input_dir, "10x_nuclei_MOp.csv")),
        .id = "dataset"
    )

    # reorder cells to match data
    row_match = match(cell_names, label_matrix$cells)
    label_matrix = label_matrix[row_match,]
    
    return(label_matrix)
}

load_single_merge_labels = function(cell_names, data_path = "../../data") {
    input_dir = file.path(data_path, "singleMerge")
    result = bind_rows(
        zeng_smart_cells = read.csv(file.path(input_dir, "SMARTer_cells_MOp_singleTree.csv")),
        zeng_smart_nuclei = read.csv(file.path(input_dir, "SMARTer_nuclei_MOp_singleTree.csv")),
        .id = "dataset"
    )
    
    # NOTE: no cell ids, we assume that the order of cells is same as data
        
    # restrict to steps where both datasets have at least 2 clusters
    n_clusters_cells = apply(result[result$dataset == "zeng_smart_cells",], 2, function(x) length(table(x)))
    n_clusters_nuclei = apply(result[result$dataset == "zeng_smart_nuclei",], 2, function(x) length(table(x)))
    keep = n_clusters_cells > 1 & n_clusters_nuclei > 1
    keep[1] = TRUE
    result = result[, keep]

    # rename columns for compatibility with analysis/visualization modules
    # (format METHOD_NAME.MERGING_STEP)
    methods = colnames(result)[-1]
    method_name = stringr::word(methods, 1, sep=stringr::fixed("."))
    method_level = stringr::word(methods, 3, sep=stringr::fixed("."))
    colnames(result)[-1] = paste(method_name, method_level, sep = ".")
    result = add_column(result, cells = cell_names, .after = 1)                          
    
    return(result)
}
                              
load_single_seurat_labels = function(cell_names, data_path = "../../data") {
    input_dir = file.path(data_path, "singleMerge")
    result = bind_rows(
        zeng_smart_cells = read_single_method(file.path(input_dir, "SMARTer_cells_MOp_singleSeurat.csv")),
        zeng_smart_nuclei = read_single_method(file.path(input_dir, "SMARTer_nuclei_MOp_singleSeurat.csv")),
        .id = "dataset"
    )

    # restrict to steps where both datasets have at least 2 clusters
    n_clusters_cells = apply(result[result$dataset == "zeng_smart_cells",], 2, function(x) length(table(x)))
    n_clusters_nuclei = apply(result[result$dataset == "zeng_smart_nuclei",], 2, function(x) length(table(x)))
    keep = n_clusters_cells > 1 & n_clusters_nuclei > 1
    keep[1] = TRUE
    result = result[, keep]

    # rename columns for compatibility with analysis/visualization modules
    # (format METHOD_NAME.PARAMETER)
    methods = colnames(result)[-(1:2)]
    method_name = "Seurat"
    method_level = gsub("\\.", "\\_", methods)
    colnames(result)[-(1:2)] = paste(method_name, method_level, sep = ".")

    # reorder cells to match data
    row_match = match(cell_names, result$cells)
    result = result[row_match,]

    return(as.data.frame(result))
}
                              
read_single_method = function(filename) {
    result = read.table(filename, check.names = FALSE)
    result = as_tibble(result, rownames = "cells")
}
