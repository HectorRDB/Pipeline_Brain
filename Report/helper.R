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

ARItrend <- function(merger) {
  clusters <- merger$initalMat
  i <- which(colnames(clusters) == "RsecT")
  clusters <- clusters[, -i]
  combo <- combn(seq_len(ncol(clusters)), 2) %>% as.data.frame()
  aris <- sapply(combo, function(x) {
    mclust::adjustedRandIndex(clusters[, x[1]], clusters[, x[2]])
  })
  merges <- merger$merges
  for (i in seq_len(nrow(merges))) {
    j <- merges[i, 1]
    clus <- clusters[, j]
    clus[clus %in% merges[i, 2:3]] <- min(merges[i, 2:3])
    clusters[, j] <- clus
    aris <- rbind(aris,
                  sapply(combo, function(x) {
      mclust::adjustedRandIndex(clusters[, x[1]], clusters[, x[2]])
      })
    )
  }
}