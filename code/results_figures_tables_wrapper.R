
# libraries 
library(haven)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyverse)

# paths
dbox <- '/Users/jordanrosenthalkay/Dropbox/Research/WLCP_replication'
output_path <- '/Users/jordanrosenthalkay/Dropbox/Research/WLCP_replication/output'
dbox <- '/Users/jordanrosenthalkay/Dropbox/Research/WLCP_replication'
path <- file.path(output_path,'model_output')
savepath <- file.path(output_path,'figures')

code_dir <- file.path(dbox, "/code/results_figures")

files <- list.files(code_dir, pattern = "\\.R$", full.names = TRUE)
for (f in sort(files)) {
  message(">> Running: ", basename(f))
  sys.source(f, envir = environment()) %>% suppressMessages()
}