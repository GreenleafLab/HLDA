#!/usr/bin/env Rscript

########################################
# Predict positive ctrl sites
########################################

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(readr)
  library(scales)
  library(rjson)
  library(ComplexHeatmap)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/plots"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = TRUE, recursive = TRUE)
setwd(wd)

# Chrombpnet models directory
modeldir <- "/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/chrombpnet_models"

# Block scientific notations
#options(scipen = 999)

# Load ArchR Project
#atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1a_atac_preprocess/lda")
#rna_obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1b_rna_preprocess/lda.rds")

##########################################################################################
# Collect all model training results
##########################################################################################

# Grab all json files for chrombpnet training metric
json.paths <- list.files(path = modeldir, pattern = "chrombpnet_metrics.json", recursive = T, full.names = T)

get_ith_element <- function(path, i) {
  elements <- unlist(strsplit(path, "/"))
  return(elements[i])
}
names <- sapply(json.paths, get_ith_element, i = 11) %>% unname()
names(json.paths) <- names

# Read the JSON files and extract metrics
require(rjson)
data.list <- list()
for (path in json.paths) {
  # Read JSON file
  json <- fromJSON(file = path)
  name <- get_ith_element(path, 11) #grab name
  split_names <- str_split(as.character(name), "\\.")
  celltype <- split_names[[1]][1]
  # Grab metrics
  df <- data.frame(
    modelname = name,
    celltype = celltype,
    counts_spearmanr = json$counts_metrics$peaks$spearmanr,
    counts_pearsonr = json$counts_metrics$peaks$pearsonr,
    counts_mse = json$counts_metrics$peaks$mse,
    profile_median_jsd = json$profile_metrics$peaks$median_jsd,
    profile_median_norm_jsd = json$profile_metrics$peaks$median_norm_jsd
    )
  data.list[[name]] <- df
}
metrics <- do.call(rbind, data.list)

write_tsv(metrics, file = paste0(wd, "/chrombpnet_model_training_metrics.tsv"))

##########################################################################################
# Plot chrombpnet training metrics
##########################################################################################
metrics <- read_tsv(paste0(wd, "/chrombpnet_model_training_metrics.tsv"))
metrics$celltype <- factor(metrics$celltype, levels = fineOrder %>% unlist() %>% unname())
covarLabel = "celltype"

pdf(paste0(wd, "/chrombpnet_model_training_metrics_boxplot.pdf"), w = 8, h = 3)
ggplot(metrics, aes(x = celltype, y = counts_spearmanr, fill = celltype)) + 
  geom_boxplot() +
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.4, 0.8) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(metrics, aes(x = celltype, y = counts_pearsonr, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.4, 0.85) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(metrics, aes(x = celltype, y = counts_mse, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.1, 0.9) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(metrics, aes(x = celltype, y = profile_median_jsd, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.2, 0.8) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(metrics, aes(x = celltype, y = profile_median_norm_jsd, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0, 0.6) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf("/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/plots/Figure02S_chrombpnet_model_training_metrics_boxplot.pdf", w = 8, h = 3)
ggplot(metrics, aes(x = celltype, y = counts_spearmanr, fill = celltype)) + 
  geom_boxplot() +
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.4, 0.8) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", panel.grid = element_blank())

ggplot(metrics, aes(x = celltype, y = counts_pearsonr, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.4, 0.85) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", panel.grid = element_blank())

ggplot(metrics, aes(x = celltype, y = counts_mse, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.1, 0.9) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", panel.grid = element_blank())

ggplot(metrics, aes(x = celltype, y = profile_median_jsd, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0.2, 0.8) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", panel.grid = element_blank())

ggplot(metrics, aes(x = celltype, y = profile_median_norm_jsd, fill = celltype)) + 
  geom_boxplot() + 
  scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name = covarLabel, na.value="grey") +
  scale_fill_manual(values=FineNamedClustCmap) +
  ylim(0, 0.6) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", panel.grid = element_blank())

dev.off()
