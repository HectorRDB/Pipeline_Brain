suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-d", "--dataset"),
              action = "store", default = "SMARTer_cells_MOp", type = "character",
              help = "Where to store the output"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$d)) {
  dataset <- opt$dataset
} else {
  stop("Missing d argument")
}

Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")


rmarkdown::render('dataset_analysis.Rmd',
                  params = list(dataset = opt$d,
                                title = paste0('Analysis of the ', opt$d,
                                               ' dataset')),
                  output_file = paste0(opt$d, '_Analysis.html'))