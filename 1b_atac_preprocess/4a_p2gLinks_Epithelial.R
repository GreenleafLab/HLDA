#!/usr/bin/env Rscript

#############################
# Peak2Gene linkages for each compartments
#############################

#Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# set working directory
subgroup <- "Epithelial"
wd <- sprintf("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/%s", subgroup)
full_dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_v2"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
addArchRGenome("hg38")
pointSize <- 1.0

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)

# Get peaks that were called on this subproject's subclusters from full ArchR project
full_proj <- loadArchRProject(full_dir, force=TRUE)
full_peaks <- getPeakSet(full_proj)
peaks <- getClusterPeaks(full_proj, clusterNames=unique(atac_proj$FineNamedClust), peakGR=full_peaks, groupBy = "FineNamedClust")
rm(full_proj); gc()

plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

###########################################################################################
# Do not proceed prior to calling peaks
###########################################################################################

# Compute group coverages (do we need to recalculate group coverages on subproject?)
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="FineNamedClust", 
  minCells = 40, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
  )

# Now add these peaks to the subproject and generate peak matrix
atac_proj <- addPeakSet(atac_proj, peakSet=peaks, force=TRUE)
atac_proj <- addPeakMatrix(atac_proj, force=TRUE)

# Load and convert PFMS from FigR paper into PWMs
cisbp_pwms <- readRDS("/oak/stanford/groups/wjg/skim/resources/MotifPWMs/FigR_cisbp_pwms/cisBP_human_pfms_2021.rds")
cisbp_pwms <- TFBSTools::toPWM(cisbp_pwms)

# Add motif annotations based on Buenrostro lab's cleaned up cisbp pwms
atac_proj <- addMotifAnnotations(atac_proj, name = "Motif", motifPWMs = cisbp_pwms, force=TRUE)

# Add Seurat PCA to ArchR Project
# obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# # Import Seurat PCA embeddings for RNA
# pca <- Embeddings(obj, reduction = "pca")
# str_sub(rownames(pca), -19, -19) <- "#" # convert cell names for ArchR
# pca <- pca[getCellNames(atac_proj),]

# # add RNA PCs to ArchR Project
# atac_proj@reducedDims[["PCA"]] <- SimpleList(
#   matDR = pca,
#   date = Sys.time(),
#   assay.used = obj@reductions$pca@assay.used,
#   scaleDims = NA,
#   corToDepth = NA
# )

# Calculate coaccessibility
atac_proj <- addCoAccessibility(
    ArchRProj = atac_proj,
    reducedDims = "PCA"
)

# Calculate peak-to-gene links
atac_proj <- addPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "PCA",
    useMatrix = "GeneExpressionMatrix"
)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

# Add UMAP coordinates based on the Seurat UMAP
# Grab umap coordinates and adjust colnames and cell barcodes to add to archr project
# Add UMAP coordinates based on the Seurat UMAP
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")
umap <- obj@reductions$umap@cell.embeddings
CBs <- rownames(umap)
str_sub(CBs, -19, -19) <- "#"
rownames(umap) <- CBs
umap <- data.frame(umap)

# Only keep cells in the atac project 
umap <- umap %>% filter(row.names(umap) %in% getCellNames(atac_proj))
colnames(umap) <- c("custom#UMAP1", "custom#UMAP2")

# Add seurat umap coordinates into ArchR
atac_proj@embeddings$customUMAP <- SimpleList(df = umap, params = list())

# Colormaps
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds")
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")
cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds")

atac_proj <- visualizeClustering(
  atac_proj, 
  pointSize=pointSize, 
  sampleCmap=sample_cmap, 
  gest_age_cmap=gest_age_cmap, 
  clusterName = "FineNamedClust", 
  embedding = "customUMAP",
  filename = "Plot-Custom-UMAP-Sample-FineNamedClust.pdf")

# Save intermediate output
saveArchRProject(atac_proj)


###########################################################################################
