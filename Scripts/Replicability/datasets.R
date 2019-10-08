
library(SingleCellExperiment)

remove_zero_genes = function(dataset) {
  return(dataset[Matrix::rowSums(assay(dataset)) > 0, ])
}

fuse_datasets = function(dataset_list) {
    result = restrict_to_common_genes(dataset_list)
    gc()
    common_col_names = Reduce(intersect, lapply(result, function(d) colnames(colData(d))))
    for (i in seq_along(result)) {
        colData(result[[i]]) = colData(result[[i]])[, common_col_names]
    }
    result = do.call(cbind, result)
    gc()
    result$study_id = rep(names(dataset_list), sapply(dataset_list, ncol))
    return(result)
}

restrict_to_common_genes = function(data_list) {
  common_genes = Reduce(intersect, lapply(data_list, rownames))
  return(lapply(data_list, function(d) d[common_genes,]))
}
                
fuse_coldata = function(dataset_list, column_name) {
    result = lapply(dataset_list, function(d) colData(d)[, column_name])
    result = do.call(rbind, result)
    result[] = lapply(result, as.character)
    return(result)
}

cell_type_average = function(dataset, cell_types) {
  cell_type_levels = levels(as.factor(cell_types))
  result = matrix(0, nrow = nrow(dataset), ncol = length(cell_type_levels))
  rownames(result) = rownames(dataset)
  colnames(result) = cell_type_levels
  for (i in seq_along(cell_type_levels)) {
    result[, i] = Matrix::rowMeans(dataset[, cell_types == cell_type_levels[i], drop = FALSE])
  }
  return(result)
}

get_deciles = function(m) {
  quantiles = seq(0, 1, length = 11)
  return(as.numeric(cut(m, quantile(m, quantiles), include.lowest = TRUE)))
}
                       
read_sparse_csv = function(filename) {
    f = file(filename, "r")
    col_names = read_header(f)
    
    row_names = c()
    row_indices = list()
    col_indices = list()
    values = list()
    row_number = 1
    line = readLines(f, 1)
    while ((length(line) > 0)) {
        data = split_line(line)
        row_names = c(row_names, data[1])
        counts = as.numeric(data[-1])
        nnz = which(counts > 0)
        row_indices[[length(row_indices)+1]] = rep(row_number, length(nnz))
        col_indices[[length(col_indices)+1]] = nnz
        values[[length(values)+1]] = counts[nnz]
        row_number = row_number + 1
        line = readLines(f, 1)
    }
    close(f)
    
    result = Matrix::sparseMatrix(i = unlist(row_indices),
                                  j = unlist(col_indices),
                                  x = unlist(values),
                                  dims = c(length(row_names), length(col_names)),
                                  dimnames = list(row_names, col_names))
    return(result)
}

read_header = function(fileptr, sep = ",") {
    line = readLines(fileptr, 1)
    return(split_line(line, sep)[-1])
}
                       
split_line = function(line, sep = ",") {
    return(strsplit(line, split = sep, fixed=TRUE)[[1]])
}
                       
reduce_to_prefix = function(labels, valid_prefixes) {
    labels = as.character(labels)
    result = labels
    for (s in valid_prefixes) {
        result[startsWith(labels, s)] = s
    }
    return(result)
}

sample_by_label = function(labels, max_size) {
    labels = as.character(labels)
    result = rep(FALSE, length(labels))
    for (l in unique(labels)) {
        is_from_label = labels == l
        if (sum(is_from_label) < max_size) {
            result[is_from_label] = TRUE
        } else {
            result[is_from_label][sample.int(sum(is_from_label), max_size)] = TRUE
        }
    }
    return(result)
}
