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

seurat_params <- function(seurat_ARI) {
  tree <- as.dendrogram(hclust(dist(seurat_ARI)))
  inds1 <- unlist(tree[[1]])
  inds2 <- unlist(tree[[2]])
  seurat_ARI1 <- seurat_ARI[inds1, inds1]
  seurat_ARI2 <- seurat_ARI[inds2, inds2]
  param1 <- names(which.max(colMeans(seurat_ARI1)))
  param2 <- names(which.max(colMeans(seurat_ARI2)))
  print(param1)
  print(param2)
  return(c(param1, param2))
}

mergeManyPairwise <- function(clusteringMatrix, nCores = 3) {
  # Turn the matrix into a numeric matrix
  clusMat <- apply(clusteringMatrix, 2, function(x) {
    x[x != "-1"] <- as.numeric(factor(x[x != "-1"]))
    x[x == "-1"] <- -1
    x <- as.integer(x)
  })
  
  # Initialize the values
  clusters <- apply(clusMat, 2, unique)
  currentMat <- clusMat
  baseARI <- apply(clusMat, 2, function(x) {
    apply(clusMat, 2, function(y) {
      mclust::adjustedRandIndex(x, y)
    })
  })
  bestARI <- baseARI
  working <- TRUE
  merges <- NULL
  ImpARI <- NULL
  
  # Try to see if any merge would increse 
  while (working) {
    # Test all pairwise clusters to merge
    # For every cluster label list
    mergeResults <- mclapply(1:ncol(currentMat), function(whClus) {
      clus <- currentMat[, whClus]
      clusternames <- clusters[[whClus]]
      clusPairs <- combn(clusternames[clusternames != -1], 2)
      
      # For every pair of labels in that list
      deltaARI <- apply(clusPairs, 2, function(pair) {
        sapply((1:ncol(clusMat))[-whClus], function(otherClus) {
          clus[clus %in% pair] <- max(clus) + 1
          mclust::adjustedRandIndex(clus, currentMat[, otherClus])
        })
      }) - bestARI[whClus, -whClus]
      
      return(colMeans(deltaARI))
    }, mc.cores = nCores)
    
    # Find best pair to merge
    maxs <- sapply(mergeResults, max)
    
    # Only merge if it improves ARI
    if (max(maxs) > 0) {
      whClus <- which.max(maxs)
      # update clusters
      clusternames <- clusters[[whClus]]
      clusPairs <- combn(clusternames[clusternames != -1], 2)
      pair <- clusPairs[, which.max(mergeResults[[whClus]])]
      indsPair <- which(currentMat[, whClus] %in% pair)
      currentMat[indsPair, whClus] <- min(pair)
      clusters[[whClus]] <- unique(currentMat[, whClus])
      
      # update bestARI
      newARIs <- sapply((1:ncol(currentMat))[-whClus], function(j) {
        mclust::adjustedRandIndex(currentMat[, whClus], currentMat[, j])
      })
      bestARI[whClus, -whClus] <- newARIs
      bestARI[-whClus, whClus] <- newARIs
      
      # tracking
      merges <- rbind(merges, c(whClus, pair))
      ImpARI <- c(ImpARI, mean(bestARI[upper.tri(bestARI)]))
      print(c(whClus, pair))
    } else {
      working <- FALSE
    }
    
    # If no more to merge in any of them, stop
    if (sum(sapply(clusters, length) == 1) == length(clusters)) stop()
  }
  
  colnames(merges) <- c("clustering", "cluster1", "cluster2")
  return(list("initalMat" = clusteringMatrix,
              "currentMat" = currentMat,
              "merges" = merges,
              "ImpARI" = ImpARI))
}
