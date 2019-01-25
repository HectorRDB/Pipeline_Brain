library(clusterExperiment)

dataset <- "SMARTer_cells_MOp/"
source("3a-zinb.R")


sequential <- FALSE
subsample <- T
clusterFunction <- "pam"

NCORES <- 
##---- load data

# load the zinbwave results

seAllen <- readRDS(file.path(dataOutput,"seAllen_zinbwave.rds"))

##---- run RSEC
zinbDim <- ncol(reducedDim(seAllen, "zinbwave10"))
print(system.time(seAllen <- RSEC(seAllen, k0s = seq(5,50,by=5), alphas = c(0.1,0.3),
                                  reduceMethod = "zinbwave10",
                                  nReducedDims = zinbDim,
                                  sequential=sequential,subsample=subsa$
                                    betas = c(0.8),
                                  clusterFunction = clusterFunction, minSizes=1,
                                  ncores = NCORES, isCount=FALSE,
                                  dendroReduce="zinbwave10", dendroNDims=zinbDim,
                                  subsampleArgs = list(resamp.num=50,
                                                       clusterFunction="kmeans"),
                                  verbose=TRUE,
                                  consensusProportion = 0.7,
                                  mergeMethod = "adjP", mergeCutoff=0.1,
                                  mergeLogFCcutoff=1, random.seed=23578$
                                    consensusMinSize = 10, run = TRUE)))

saveRDS(seAllen, file = file.path(dataOutput,sprintf("seAllen_RSEC_%s.rds",tag)))