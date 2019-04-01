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

plotARIReduce <- function(merger) {
  # Before, No unclustered cells for RSEC
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  InitialARI <- apply(merger$initalMat[, -r1], 2, function(x) {
    apply(merger$initalMat[, -r1], 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  
  p1 <- plotARIs(InitialARI) +
    ggtitle("ARI before any merging, all cells assigned") +
    theme(title = element_text(size = 6))
  
  # Before, All cells assigned
  r2 <- which(colnames(merger$initalMat) == "Rsec")
  InitialARI <- apply(merger$initalMat[,-r2], 2, function(x) {
    apply(merger$initalMat[,-r2], 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  
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

assignRsec <- function(merger) {
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  r2 <- which(colnames(merger$initalMat) == "Rsec")
  
  currentMat <- merger$currentMat
  Rsec_merges <- merger$merges
  Rsec_merges <- Rsec_merges[Rsec_merges[,1] == 2, ]
  Rsec_merges <- Rsec_merges[, -1]
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

type <- function(dataset) {
  if (str_detect(dataset, "SMART")) return("Smart-Seq")
  if (str_detect(dataset, "10x")) return("10X")
  stop("Type unknown")
}

#' ARI improvement
#' 
#' Compute the ARI improvement over the ARI merging procedure
#' @param merger the result from having run mergeManyPairwise on the dataset
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
#' @param merger the result from having run mergeManyPairwise on the dataset
#' @return a ggplot object
ARItrend <- function(merger) {
 baseMat <- merger$initalMat
 j <- which(colnames(baseMat) == "RsecT")
 baseMat <- baseMat[, -j]
 # baseARI <- apply(baseMat, 2, function(x) {
 #   apply(baseMat, 2, function(y) {
 #     adjustedRandIndex(x, y)
 #   })
 # })
 # baseARI <- baseARI[upper.tri(baseARI)] %>% mean()
 # ARI <- c(baseARI, merger$ImpARI)
 ARI <- ARIImp(merger)
 n_clus <- lapply(1:nrow(merger$merges), function(m){
              diff <- rep(0, 3)
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
   mutate(type = ifelse(change == "ARI_Imp", "ARI Improvement", "Cluster size"))
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
  p
}

#' Find the clustering matrix that we would get if we stopped the ARI merging 
#' early
#' @param  merger the result from having run mergeManyPairwise on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement
#' @return A matrix with the same dimensions as the currentmMat of the merger
#' argument
intermediateMat <- function(merger, p = .9) {
  # Compute ARI imp and find where to stop the merge
  ARI <- ARIImp(merger)
  int_merges <- merger$merges
  j <- min(which(ARI >= min(ARI) + p * (max(ARI) - min(ARI))))
  int_merges <- int_merges[1:j, ]
  assign <- sapply(1:ncol(merger$currentMat), function(clus) {
    if (sum(int_merges[, 1] == clus) == 0) {
      return(merger$initalMat[,clus])
    } else {
      clus_merges <- int_merges[int_merges[, 1] == clus, ] %>%
        as.matrix(matrix(ncol = 3))
      lapply(1:nrow(merger$currentMat), function(i){
        cell <- merger$initalMat[i, clus]
        for (j in 1:nrow(clus_merges)) {
          if (cell %in% clus_merges[j, ]) cell <- min(clus_merges[j, ])
        }
        return(cell)
      }) %>% unlist() %>% return()
    }
  }) 
  
  j <- which(colnames(merger$initalMat) == "RsecT")
  colnames(assign) <- colnames(merger$initalMat)[-j]
  return(assign)
}
