suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the object after running"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc, "\n")
} else {
  stop("Missing l argument\n")
}
if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument\n")
}

read_table(loc) %>%
  filter(X1 == "Mem:") %>%
  mutate(used = str_remove_all(used, "G") %>% as.numeric()) %>%
  summarise(Max = max(used),
            mean = mean(used)) %>%
  write.table(output)
