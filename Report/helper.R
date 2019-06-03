#' Plot an heatmap of the ARI matrix
#' 
#' We can compute the ARI between pairs of cluster labels. This function plots
#' a matrix where a cell is the adjusted Rand Index between cluster label of
#' row i and cluster label of column j. 
#' @param ARI the matrix of pairwise ARI
#' @param small whether to also plot the values
#' @return a ggplot object

plotARIs <- function(ARI, small = T) {
  p <- ARI %>% as.data.frame() %>%
    mutate(label = rownames(ARI)) %>%
    gather(key = label2, value = ari, -(ncol(ARI) + 1)) %>%
    ggplot(aes(x = label, y = label2, fill = ari)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme_classic() +
    theme(axis.line = element_blank())
  if (small) {
    p <- p  +
      geom_text(aes(label = round(ari, 2))) +
      guides(fill = F)
  }
  return(p)
}

#' Plot the reduction in cluster size for an ARI merging
#' @param merger The output from an ARI merging
#' @return a ggplot object
plotPrePost <- function(merger) {
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  pre <- apply(merger$initalMat[,-r1], 2, function(x) length(unique(x)))
  post <- apply(merger$currentMat, 2, function(x) length(unique(x)))
  df <- data.frame(methods = names(pre),
                   before = pre,
                   after = post) %>%
    gather(key = "time", value = "Nb", -methods) %>%
    mutate(time = factor(time, levels = c("before", "after")))
  p <- ggplot(df, aes(x = methods, y = Nb, fill = time)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    scale_fill_viridis_d(option = "E") +
    labs(x = "Clustering Methods", y = "Number of clusters", fill = "") +
    ggtitle("Reduction in number of clusters with ARI merging")
  return(p)
}

#' plot the ARI improvement between methods
#' 
#' The output from this function is a grid of 4 plots, based on the
#' \code{\link{plotARIs}} function. The first row of 2 plots is before the 
#' merging procedure, the second is after the merging procedure. The first 
#' column is when RSEC uses all assigned cells, the second column is when only
#' using the cells that RSEC cluster for computing the ARI between RSEC and 
#' another partition of the data.
#' @param merger the result from having run \code{\link{mergeManyPairwise}} 
#' on the dataset
#' @return the output from \code{\link{plot_grid}}
#' @import cowplot
plotARIReduce <- function(merger) {
  # Before, No unclustered cells for RSEC
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  InitialARI <- apply(merger$initalMat[, -r1], 2, function(x) {
    apply(merger$initalMat[, -r1], 2, function(y) {
      inds <- x != -1 & y != -1
      xa <- x[inds]
      ya <- y[inds]
      adjustedRandIndex(xa, ya)
    })
  })
  
  p1 <- plotARIs(InitialARI) +
    ggtitle("ARI before any merging, no unclustered cells for RSEC") +
    theme(title = element_text(size = 6))
  
  # Before, All cells assigned
  r2 <- which(colnames(merger$initalMat) == "Rsec")
  InitialARI <- apply(merger$initalMat[,-r2], 2, function(x) {
    apply(merger$initalMat[,-r2], 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  colnames(InitialARI)[colnames(InitialARI) == "RsecT"] <- "Rsec"
  rownames(InitialARI)[rownames(InitialARI) == "RsecT"] <- "Rsec"
  
  p2 <- plotARIs(InitialARI) +
    ggtitle("ARI before any merging, all cells assigned") +
    theme(title = element_text(size = 6))
  
  ## After, No unclustered cells for Rsec
  FinalARI <- apply(merger$currentMat, 2, function(x) {
    apply(merger$currentMat, 2, function(y) {
      inds <- x != -1 & y != -1
      xa <- x[inds]
      ya <- y[inds]
      adjustedRandIndex(xa, ya)
    })
  })
  
  p3 <- plotARIs(FinalARI) +
    ggtitle("ARI after merging, no unclustered cells for RSEC") +
    theme(title = element_text(size = 6))
  
  ## After, Rsec with all cells
  currentMat <- merger$currentMat
  currentMat[, "Rsec"] <- assignRsec(merger) 
  FinalARI <- apply(currentMat, 2, function(x) {
    apply(currentMat, 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  
  p4 <- plotARIs(FinalARI) +
    ggtitle("ARI after merging, all cells assigned") +
    theme(title = element_text(size = 6))
  
  plot_grid(ggdraw() + draw_plot(p1),
            ggdraw() + draw_plot(p2),
            ggdraw() + draw_plot(p3),
            ggdraw() + draw_plot(p4),
            ncol = 2, rel_heights = rep(.25, 4), rel_widths = rep(.25, 4))
}

#' Assign cells using the assignUnassigned function of RSEC
#' 
#' By default, RSEC does not assign all cells, leaving those as "-1". 
#' However, to give a fair comparison with other labels, it is necessary to 
#' assign cells so it is necessary to track that along the merging
#' @param merger the result from having run \code{\link{mergeManyPairwise}} 
#' on the dataset
#' @param p when to stop the merging, when mean ARI has improved to p (between 0
#' and 1) of the final value.
assignRsec <- function(merger, p = 1) {
  ARI <- ARIImp(merger)
  K <- min(which(ARI >= min(ARI) + p * (max(ARI) - min(ARI))))
  
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  r2 <- which(colnames(merger$initalMat) == "Rsec")
  
  currentMat <- merger$currentMat
  Rsec_merges <- merger$merges[1:(K - 1), ]
  Rsec_merges <- Rsec_merges[Rsec_merges[,1] == 2, ]
  if (is.null(dim(Rsec_merges))) Rsec_merges <- matrix(Rsec_merges, nrow = 1)
  Rsec_merges <- Rsec_merges[, -1]
  if (is.null(dim(Rsec_merges))) Rsec_merges <- matrix(Rsec_merges, nrow = 1)
  if (nrow(Rsec_merges) == 0) {
    return(merger$initalMat$RsecT)
  } else {
    assign <- lapply(1:nrow(currentMat), function(i) {
      cell <- merger$initalMat[i, r2]
      cellT <- merger$initalMat[i, r1]
      for (j in 1:nrow(Rsec_merges)) {
        if (cellT %in% Rsec_merges[j, ]) cellT <- min(Rsec_merges[j, ])
      }
      return(cellT)
    }) %>%
      unlist()
    return(assign)
  }
}

type <- function(dataset) {
  if (str_detect(dataset, "SMART")) return("Smart-Seq")
  if (str_detect(dataset, "10x")) return("10X")
  stop("Type unknown")
}

#' ARI improvement
#' 
#' Compute the ARI improvement over the ARI merging procedure
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#'  on the dataset
#' @return a vector with the mean ARI between methods at each step
ARIImp <- function(merger) {
  baseMat <- merger$initalMat
  j <- which(colnames(baseMat) == "RsecT")
  baseMat <- baseMat[, -j]
  baseARI <- apply(baseMat, 2, function(x) {
    apply(baseMat, 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  baseARI <- baseARI[upper.tri(baseARI)] %>% mean()
  ARI <- c(baseARI, merger$ImpARI)
  ARI <- cumsum(ARI)
  return(ARI)
}

#' ARI improvement plot
#' 
#' A plot to see how ARI improves over merging
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#'  on the dataset
#' @return a ggplot object
ARItrend <- function(merger) {
 baseMat <- merger$initalMat
 j <- which(colnames(baseMat) == "RsecT")
 baseMat <- baseMat[, -j]
 ARI <- ARIImp(merger)
 n_clus <- lapply(1:nrow(merger$merges), function(m){
              diff <- rep(0, ncol(baseMat))
              diff[merger$merges[m, 1]] <- -1
              matrix(diff, nrow = 1)
             }) %>%
   do.call('rbind', args = .)
 n_clus <- rbind(sapply(baseMat, n_distinct) %>% matrix(data = ., nrow = 1),
                 n_clus)
 n_clus <- apply(n_clus, 2, cumsum)
 colnames(n_clus) <- colnames(baseMat)
 df <- data.frame(step = 0:length(merger$ImpARI),
                  ARI_Imp = ARI,
                  n_clus) %>%
   gather(key = "change", value = "value", -step) %>%
   mutate(type = ifelse(change == "ARI_Imp", "ARI Improvement", "Nb of clusters"))
 p <- ggplot(df, aes(x = step, y = value)) +
   geom_path(size = 2, aes(group = change, col = change)) +
   facet_wrap(~type, scales = "free") +
   theme_classic() +
   scale_x_continuous(breaks = c(0, length(merger$ImpARI)),
                      labels = c("Initial", "Final")) +
   geom_hline(yintercept = min(ARI) + .9 * (max(ARI) - min(ARI)),
              col = "grey", linetype = "dashed", size = 2) +
   geom_vline(xintercept = min(which(ARI >= min(ARI) +
                                       .9 * (max(ARI) - min(ARI)))),
              col = "grey", linetype = "dashed", size = 2) +
   labs(y = "Change over merging",
        col = "type")
  return(p)
}

#' Find the clustering matrix that we would get if we stopped the ARI merging 
#' early
#' @param merger the result from having run \code{\link{mergeManyPairwise}} 
#' on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement
#' @return A matrix with the same dimensions as the currentmMat of the merger
#' argument
intermediateMat <- function(merger, p = .9) {
  # Compute ARI imp and find where to stop the merge
  ARI <- ARIImp(merger)
  int_merges <- merger$merges
  j <- min(which(ARI[2:length(ARI)] >= min(ARI) + p * (max(ARI) - min(ARI))))
  int_merges <- int_merges[1:j, ]
  assign <- sapply(colnames(merger$currentMat), function(clus) {
    J <- which(colnames(merger$initalMat) == clus)
    if (sum(int_merges[, 1] == J) == 0) {
      return(merger$initalMat[, clus])
    } else {
      clus_merges <- int_merges[int_merges[, 1] == J, ] %>%
        as.matrix() %>% matrix(ncol = 3)
      cells <- merger$initalMat[, clus]
      walk(1:nrow(clus_merges), function(i){
        indsPair <- which(cells %in% clus_merges[i, 2:3])
        cells[indsPair] <- min(clus_merges[i, 2:3])
      })
      return(cells)
    }
  }) 
  
  j <- which(colnames(merger$initalMat) == "RsecT")
  colnames(assign) <- colnames(merger$initalMat)[-j]
  return(assign)
}

#' Track the evolution of a function along merging 
#' 
#' For a given ARI merging, compute the evolution on the function f with
#' another partition
#' @param merger the result from having run \code{\link{mergeManyPairwise}} 
#' on the dataset
#' @param f the function used, can be computed on any partition of the space 
#' @param ... additional arguments passed to f
#' @return a matrix with a column per initial clustering, and a row per merge
#' with the f value computed
FTracking <- function(merger, f, ...){
  # Go over the merge and compute the homogeneity as we go
  baseMat <- merger$initalMat
  j <- which(colnames(baseMat) == "Rsec")
  baseMat <- baseMat[, -j]
  currentMat <- baseMat
  
  Evolution <- apply(baseMat, 2, f, ...) %>%  matrix(ncol = ncol(baseMat))
  
  for (m in seq_len(nrow(merger$merges))) {
    wClus <- merger$merges[m, 1]
    clus <- currentMat[, wClus]
    pair <- merger$merges[m, 2:3]
    clus[clus %in% pair] <- min(pair)
    currentMat[, wClus] <- clus
    Evolution <- rbind(Evolution, Evolution[nrow(Evolution), ])
    Evolution[nrow(Evolution), wClus] <- f(clus, ...)
  }
  return(Evolution)
}


#' Track the evolution of the ARI between allen cluster and the consensus
#' 
#' For a given ARI merging, compute the evolution of the ARI between the allen
#' cluster and the consensus merger.
#' @param merger the result from having run \code{\link{mergeManyPairwise}} 
#' on the dataset
#' @param allen1 the first set of allen clustering
#' @param allen2 the second set of allen clustering
#' @param verbose Whether to print a progress bar
#' @param ... Other arguments passed to \code{\link{MakeConsensus}}
#' @import progress
#' @return a vector with the ARI
ARItrendAllen <- function(merger, allen1, allen2, verbose = TRUE, ...){
  # Go over the merge and compute the homogeneity as we go
  baseMat <- merger$initalMat
  j <- which(colnames(baseMat) == "RsecT")
  baseMat <- baseMat[, -j]
  if (verbose) {
    pb <- progress_estimated(n = nrow(merger$merges) + 1)
  }
  baseConsensus <- Consensus(baseMat, ...)
  if (verbose) pb$tick()
  inds <- baseConsensus == -1
  currentMat <- baseMat
  
  ARI <- matrix(c(adjustedRandIndex(allen1[!inds], baseConsensus[!inds]),
                  adjustedRandIndex(allen2[!inds], baseConsensus[!inds])), 
                ncol = 2)
  for (m in seq_len(nrow(merger$merges))) {
    wClus <- merger$merges[m, 1]
    clus <- currentMat[, wClus]
    pair <- merger$merges[m, 2:3]
    clus[clus %in% pair] <- min(pair)
    currentMat[, wClus] <- clus
    currentConsensus <- Consensus(currentMat, ...)
    if (verbose) pb$tick()$print()
    inds <- currentConsensus == -1
    ARI <- rbind(ARI, matrix(
                  c(adjustedRandIndex(allen1[!inds], currentConsensus[!inds]),
                    adjustedRandIndex(allen2[!inds], currentConsensus[!inds])),
                  ncol = 2)
                 )
  }
  return(ARI)
}

#' Find the consensus clustering between three methods and return the consensus
#' @param clusMat The clustering matrix with a row per cell and a column per 
#' clustering label type
#' @param ... Other arguments passed to \code{\link{MakeConsensus}}
Consensus <- function(clusMat, ...) {
  cellsConsensus <- suppressWarnings(
    makeConsensus(x = as.matrix(clusMat), clusterLabel = "makeConsensus",
                  proportion = 2/3, minSize = 100, ...)
  )
  return(cellsConsensus$clustering)
}
