suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc, "\n")
} else {
  stop("Missing l argument\n")
}

usage <- read_table(loc) %>%
  filter(X1 == "Mem:") %>%
  mutate(used = str_remove_all(used, "G") %>% as.numeric())
usage %>% summarise(Max = max(used), mean = mean(used))

plot(1:nrow(usage), usage$used, type = "l")
