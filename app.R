# - Load all required packages -
# -- Bioconductor packages --
# source("https://bioconductor.org/biocLite.R")
# first time these bioconductor packages have to be installed with biocLite("phyloseq")
library(phyloseq)
library(vegan)
library(DESeq2)
# ----
library(ggplot2)
library(shiny)
library(lubridate)
library(tidyr)
library(dplyr)
library(viridis)
# library(rdrop2)
#library(xtable)
library(pheatmap)
# library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
# --



# - load functions -
functionpath <- "./Functions"
# functionpath <- "/Users/jvb740/Coursera_MOOC/20161202_LearningShiny_FantasySports/shinyy/Apps/Teaching_Apps/MicrobiomeX4_App/Functions"
source(file.path(functionpath, "_n_000_helper_functions.R"))
source(file.path(functionpath, "_n_010_explore_ps_functions.R"))
source(file.path(functionpath, "_n_020_alpha_diversity_functions.R"))
source(file.path(functionpath, "_n_030_preprocess_filtering_functions.R"))
source(file.path(functionpath, "_n_040_beta_diversity_functions.R"))
source(file.path(functionpath, "_n_050_diff_abundance_functions.R"))
source(file.path(functionpath, "_n_060_phylum_analysis_functions.R"))
source(file.path(functionpath, "shiny_app_functions.R"))
# --



# - load the ui files -
source("MicrobiomeX4_ui.R")
# --



# - load the server function -
source("MicrobiomeX4_server.R")
# --



# - run the app -
shinyApp(ui = ui, server = server)
# --