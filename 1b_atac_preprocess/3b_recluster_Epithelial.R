#!/usr/bin/env Rscript

#############################
# Subgroup Clustering
#############################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# set working directory
subgroup <- "Epithelial"
wd <- sprintf("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/%s", subgroup)
maxClusters <- c(3, 5, 7)

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

addArchRGenome("hg38")

pointSize <- 0.5

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_obj <- readRDS(paste0("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/", subgroup, "/", subgroup, ".rds"))
  
plotDir <- paste0(wd, "/Plots")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# Colormaps
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds")
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")
cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds")

##########################################################################################
# Re-cluster subclustered ArchR project
##########################################################################################

atac_proj <- addIterativeLSI(
    ArchRProj = atac_proj,
    useMatrix = "TileMatrix", 
    iterations = 4,
    name = "LSI_ATAC",
    clusterParams = list(
      resolution = c(0.1, 0.2, 0.4),
      maxClusters = maxClusters
    ),
    sampleCellsPre = 10000,
    dimsToUse = 1:15,
    force = TRUE
)

# (Batch correct sample effects)
atac_proj <- addHarmony(
    ArchRProj = atac_proj,
    reducedDims = "LSI_ATAC",
    name = "Harmony_ATAC",
    groupBy = "age",
    force = TRUE
)

# Add umaps
atac_proj <- addUMAP(
    ArchRProj = atac_proj, 
    reducedDims = "LSI_ATAC", 
    name = "UMAP_ATAC", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

atac_proj <- addUMAP(
    ArchRProj = atac_proj, 
    reducedDims = "Harmony_ATAC", 
    name = "UMAP_Harmony_ATAC", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

# Call clusters
atac_proj <- addClusters(atac_proj, reducedDims = "Harmony_ATAC", name = "Clusters_ATAC", resolution = 0.2, force = TRUE)

# Make various cluster plots:

# FineNamedClust on all the different UMAPs
atac_proj <- visualizeClustering(
  atac_proj, 
  pointSize=pointSize, 
  sampleCmap=sample_cmap, 
  gest_age_cmap=gest_age_cmap, 
  clusterName = "Clusters_ATAC", 
  embedding = "UMAP_ATAC",
  filename = "Plot-UMAP-ATAC-Sample-Clusters_ATAC.pdf")

# FineNamedClust on all the different UMAPs with harmony
atac_proj <- visualizeClustering(
  atac_proj, 
  pointSize=pointSize, 
  sampleCmap=sample_cmap, 
  gest_age_cmap=gest_age_cmap, 
  clusterName = "Clusters_ATAC", 
  embedding = "UMAP_Harmony_ATAC",
  filename = "Plot-UMAP-Harmony-ATAC-Sample-Clusters_ATAC.pdf")

saveArchRProject(atac_proj)

# Plot ATAC cluster on RNA UMAP to check for doublets or low quality cells
clusters <- unique(atac_proj$Clusters_ATAC) %>% sort()

pdf(paste0(plotDir, "/Plot-Seurat-UMAP-RNA-Clusters_ATAC.pdf"))
for (cluster in clusters) {
  cells <- atac_proj$cellNames[which(atac_proj$Clusters_ATAC %in% cluster)]
  stringr::str_sub(cells, start = -19, end = -19) <- "_"
  print(DimPlot(rna_obj, cells.highlight = cells) + ggtitle(cluster))
}
dev.off()

###########################################################################################