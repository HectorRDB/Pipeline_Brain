clusterMatToAri <- function(seurat) {
  ARI <- apply(seurat, 2, function(x) {
    apply(seurat, 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  return(plotARIs(ARI, small = TRUE))
}
