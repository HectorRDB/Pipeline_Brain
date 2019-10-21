library(here)
library(rmarkdown)
for (dataset in c("SMARTer_cells_MOp", "SMARTer_nuclei_MOp")) {
  rmarkdown::render(here("Report", 'ARI_comp.Rmd'),
                    params = list(dataset = dataset,
                                  title = paste0('Analysis of the ', dataset,
                                                 ' dataset')),
                    output_file = paste0(dataset, '_ARI_comp.html'))
}
