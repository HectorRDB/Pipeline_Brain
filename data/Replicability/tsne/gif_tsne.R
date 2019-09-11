# Inspired by https://ryanpeek.github.io/2016-10-19-animated-gif_maps_in_R/
library(magick)
library(purrr)
library(stringr)
library(dplyr)
library(here)
types <- c("zeng_10x_cells", "zeng_10x_nuclei", "zeng_smart_cells",
           "zeng_smart_nuclei")
number <- c("_Initial.png", "_33.png", "_66.png", "_90.png", "_Final.png")
walk(types, function(type){
  files <- paste0(type, number)
  map(files, function(file){
    image_read(here("Replicability", "tsne", file))
  }) %>%
    image_join() %>% # joins image
    image_animate(fps = 1) %>% # animates, can opt for number of loops
    image_write(here("Replicability", "tsne", paste0(type, ".gif")))
})
