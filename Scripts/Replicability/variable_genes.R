
variable_genes = function(dataset, sample_size = 50000, i = 1) {
  if (sample_size < ncol(dataset)) {
    subset = dataset[,sample.int(ncol(dataset), sample_size)]
  } else {
    subset = dataset
  }
  result = mn_variable_genes_legacy(subset, i, exp_labels=subset$study_id)
  return(Reduce(intersect, result))
}

variable_genes_by_dataset = function(dataset, sample_size = 50000, i = 1) {
  if (sample_size < ncol(dataset)) {
    subset = dataset[,sample.int(ncol(dataset), sample_size)]
  } else {
    subset = dataset
  }
  result = mn_variable_genes_legacy(subset, i, exp_labels=subset$study_id)
  return(result)
}

mn_variable_genes_legacy = function(dat, i = 1, exp_labels) {

  dat = SummarizedExperiment::assay(dat, i = i)
  var_genes1 = vector("list")

  #check length of exp_labels equal # of samples
  if(length(exp_labels) != length(colnames(dat))){
      stop('experiment_labels length does not match number of samples')
  }

  #check obj contains more than 1 unique study_id
  if(length(unique(exp_labels)) < 2){
      stop('Found only 1 unique exp_labels. Please use data from more than 1 study!')
  }

  experiments = unique(exp_labels)
  j = 1
  for(exp in experiments){
    var_genes1[[j]] = mn_variable_genes_one_dataset(as.matrix(dat[ , exp_labels == exp]))
    j = j+1
  }
  return(var_genes1)
}

mn_variable_genes_one_dataset = function(data_subset) {
    genes_list = vector("list")
    median_data = matrixStats::rowMedians(data_subset)
    variance_data = matrixStats::rowVars(data_subset)
    names(variance_data) = rownames(data_subset)
    quant_med = unique(
      stats::quantile(median_data, probs = seq(0, 1, length = 11), type = 5)
    )
    genes_list = vector("list", length = length(quant_med))

    for(i in seq_along(quant_med)){
      if(i == 1){
        filt1 = median_data <= quant_med[i]
      } else {
        filt1 = median_data <= quant_med[i] & median_data > quant_med[i-1]
      }
      var_temp = variance_data[filt1]
      quant_var = stats::quantile(var_temp, na.rm = TRUE)
      filt2 = var_temp > quant_var[4]
      genes_list[[i]] = names(var_temp)[filt2]
    }
    temp = length(genes_list)
    return(unlist(genes_list[1:temp-1]))
}


scran_variable_genes = function(dataset, n = 200) {
  if (length(unique(dataset$study_id)) > 1) {
    design = model.matrix(~as.character(dataset$study_id))
    fit = scran::trendVar(dataset, design=design, use.spikes=FALSE)
  } else {
    fit = scran::trendVar(dataset, use.spikes=FALSE)
  }
  result = scran::decomposeVar(dataset, fit)
  result = result[order(result$bio, decreasing=TRUE), ]
  return(rownames(result)[seq_len(n)])
}

variable_genes_lcv = function(dataset, window_half_size = 50, number_datasets = NULL) {
    result = list()
    for (exp in unique(dataset$study_id)) {
        data_subset = as.matrix(assay(dataset[, dataset$study_id == exp]))
        
        gene_detection = rowSums(data_subset > 0)
        quantiles = quantile(gene_detection[gene_detection > 0], probs = c(0.1,0.9))
        data_subset = data_subset[gene_detection > quantiles[1] & gene_detection < quantiles[2], ]
        
        gene_order = order(rowMeans(log1p(data_subset)))
        cv_data = matrixStats::rowSds(data_subset) / rowMeans(data_subset)
        cv_data = cv_data[gene_order]
        
        var_genes = rep(0, length(cv_data))
        names(var_genes) = names(cv_data)
        for (i in seq_along(cv_data)) {
            m = max(1, i - window_half_size)
            M = min(length(cv_data), i + window_half_size)
            r = rank(cv_data[m:M], ties.method = "average") / (M-m+1)
            var_genes[i] = r[i-m+1]
        }
        var_genes = sort(var_genes, decreasing = TRUE)
        result[[length(result)+1]] = names(var_genes[seq_len(length(var_genes)*0.25)])
    }
    result = table(unlist(result))
    if (is.null(number_datasets)) {
        number_datasets = max(result)
    }
    result = names(result)[result >= number_datasets]
    return(result)
}
