suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-d", "--dataset"),
              action = "store", default = "SMARTer_cells_MOp", type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-l", "--local"),
              action = "store", default = "TRUE", type = "logical",
              help = "Whether we are running the script locally or on xsede"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$d)) {
  dataset <- opt$dataset
} else {
  stop("Missing d argument")
}

if (opt$l) {
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")  
} else {
  .libPaths("/pylon5/ib5phhp/shared/rpack/3.5")
}

rmarkdown::render('dataset_analysis.Rmd',
                  params = list(dataset = dataset,
                                title = paste0('Analysis of the ', dataset,
                                               ' dataset')),
                  output_file = paste0(dataset, '_Analysis.html'))
  