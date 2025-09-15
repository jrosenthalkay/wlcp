
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
library(countrycode)

# --- paths -----
# >>>>>>>> SET YOUR ROOT FOLDER HERE <<<<<<<<
root <- 'C:/Users/path/to/root/wlcp'
# ^ Change this path to the root of the replication package on your machine

output_path <- file.path(root,'output')
path <- file.path(root,'data/model_output')
savepath <- file.path(output_path,'figures')
code_dir <- file.path(root, "/code/results_figures")

# --- loop over all figures files -----
files <- list.files(code_dir, pattern = "\\.R$", full.names = TRUE)
for (f in sort(files)) {
  message(">> Running: ", basename(f))
  sys.source(f, envir = environment()) %>% suppressMessages()
}