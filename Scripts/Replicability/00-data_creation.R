suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(Seurat)
})

source("variable_genes.R")
# source("identifier_conversion.R")
source("datasets.R")

create_data = function() {
    dataset = list(
        smart_cells = readRDS("~/data/biccn/parsed_data/zeng_smart_cells.rds"),
        smart_nuclei = readRDS("~/data/biccn/parsed_data/zeng_smart_nuclei.rds"),
        tenx_v2_cells = readRDS("~/data/biccn/parsed_data/zeng_10x_cells.rds"),
        tenx_v2_nuclei = readRDS("~/data/biccn/parsed_data/zeng_10x_nuclei.rds")
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
    saveRDS(dataset, here("data", "full_data.rds"))
}

create_lab_data <- function() {
    counts_Zeng <- readRDS(here("data", "full_data.rds"))
    counts_Zeng <- counts_Zeng[, counts_Zeng$study_id == "zeng_10x_nuclei"]
    counts_Zeng <- assays(counts_Zeng)[[1]]
    counts_Zeng[is.na(counts_Zeng)] <- 0
    counts_Zeng <- as.matrix(counts_Zeng)
    
    counts_Regev <- read.csv("/pylon5/ib5phhp/hectorrb/Regev/count_matrix.csv",
                             row.names = 1)
    genes <- intersect(rownames(counts_Regev), rownames(counts_Zeng))
    print(sum(counts_Regev[!rownames(counts_Regev) %in% genes, ])/sum(counts_Regev))
    counts <- cbind(counts_Zeng[genes, ], counts_Regev[genes, ])
    meta <- data.frame(cells = colnames(counts))
    meta$study_id <- c(rep("Zeng", ncol(counts_Zeng)),
                       rep("Regev", ncol(counts_Regev)))
    dataset <- SingleCellExperiment(assays = list(as.matrix(counts)),
                                    colData = meta)
    
    hvg <- variable_genes(dataset)
    hvg <- intersect(rownames(dataset), hvg)
    
    rowData(dataset)$is_hvg <- rownames(dataset) %in% hvg
    dataset$class_label <- "1"
    saveRDS(dataset, here("data", "lab_data.rds"))
}

create_lab_data()