#!/usr/bin/env Rscript

#####################################
# Create color map objects
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(BuenColors)
})

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/final"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Sample color map
sample_cmap <- c("#99b8ba", "#839dad", "#7890a7", "#6d82a0", "#62759a", "#4c5a8d", "#414d86", "#2b3279", "#202573", "#0a0a66")
names(sample_cmap) <- c("PCW12", "PCW14", "PCW15.D", "PCW15.P", "PCW16", "PCW19.LU", "PCW19.RL", "PCW21.D", "PCW21.P", "PCW23")
saveRDS(sample_cmap, paste0(wd, "/sample_cmap.rds"))

#Gestational age color map
gest_age_cmap <- c("#FFF7BC", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506")
names(gest_age_cmap) <- c("PCW12", "PCW14", "PCW15", "PCW16", "PCW19", "PCW21", "PCW23")
saveRDS(gest_age_cmap, paste0(wd, "/gest_age_cmap.rds"))
