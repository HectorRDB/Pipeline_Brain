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
    image_read(here("Figures", "tSNE", file))
  }) %>%
    image_join() %>% # joins image
    image_animate(fps = 1) %>% # animates, can opt for number of loops
    image_write(here("Figures", "tSNE", paste0(type, ".gif")))
})

types <- c("10x_cells", "10x_nuclei", "SMARTer_cells", "SMARTer_nuclei")
walk(types, function(type){
  files <- list.files(here("Figures", "EDA")) %>%
    str_subset(pattern = type) %>%
    str_subset(pattern = "gif", negate = TRUE)
  map(files, function(file){
    K <- str_extract(file, "..\\.png$") %>% str_remove("\\.png$")
    image_read(here("Figures", "EDA", file)) %>%
      image_annotate(., text = paste0("K = ", K), gravity = "north", size = 50) %>%
      return()
  }) %>%
    image_join() %>% # joins image
    image_animate(fps = 1) %>% # animates, can opt for number of loops
    image_write(here("Figures", "EDA", paste0(type, "_MOp.gif")))
})
