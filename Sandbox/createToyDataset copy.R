library(here)
library(tidyverse)
# Tool dataset

set.seed(24)
init <- sample(1:10, 100, replace = TRUE)
clusMat <- matrix(rep(init, 5), ncol = 5, byrow = FALSE)
clusMat[clusMat[, 2] == 8, 2] <- sample(11:12, sum(clusMat[, 2] == 8), TRUE)
clusMat[clusMat[, 3] == 8, 3] <- sample(11:12, sum(clusMat[, 3] == 8), TRUE)
clusMat[clusMat[, 2] == 4, 4] <- sample(11:13, sum(clusMat[, 4] == 4), TRUE)
clusMat[clusMat[, 2] == 2, 5] <- sample(11:14, sum(clusMat[, 5] == 2), TRUE)