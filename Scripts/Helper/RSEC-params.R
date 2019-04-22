libs <- c("here", "tidyverse")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
load(here("Scripts", "Helper", "RSEC-params.rda"))
df1 <- params1$paramMatrix
df2 <- params2$paramMatrix
