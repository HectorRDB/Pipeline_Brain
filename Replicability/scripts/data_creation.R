suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(tidyverse)
})

current_dir = getwd()
common_dir = "~/projects/common/"

setwd(common_dir)
source("variable_genes.R")
source("identifier_conversion.R")
source("datasets.R")
setwd(current_dir)


create_data = function() {
    dataset = list(
        smart_cells = readRDS("~/data/biccn/parsed_data/zeng_smart_cells.rds"),
        smart_nuclei = readRDS("~/data/biccn/parsed_data/zeng_smart_nuclei.rds"),
        tenx_cells = readRDS("~/data/biccn/parsed_data/zeng_10x_cells.rds"),
        tenx_nuclei = readRDS("~/data/biccn/parsed_data/zeng_10x_nuclei.rds")
    )

    rownames(dataset$tenx_cells) =  convert_to_mgi_symbols_from_10x(rownames(dataset$tenx_cells))
    rownames(dataset$tenx_nuclei) = convert_to_mgi_symbols_from_10x(rownames(dataset$tenx_nuclei))

    col_data = fuse_coldata(dataset, c("cluster_label", "subclass_label", "class_label", "study_id"))
    dataset = fuse_datasets(dataset)
    colData(dataset) = col_data

    hvg = readRDS("~/projects/biccn/results/allen_broad/hvg.rds")
    hvg = convert_to_mgi_symbols_from_10x(hvg)
    hvg = intersect(rownames(dataset), hvg)
    
    rowData(dataset)$is_hvg = rownames(dataset) %in% hvg
    saveRDS(dataset, "full_data.rds")
}