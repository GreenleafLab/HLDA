#!/usr/bin/env Rscript

#####################################################################
# Call peaks on clusters from overall clustering
#####################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(Seurat)
  library(tidyr)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess"

# New directory
outdir <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_v2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

#Load Genome Annotations
addArchRGenome("hg38")
pointSize <- 0.25

##########################################################################################
# Load Previously Prepared ArchR project
##########################################################################################

orig_proj <- loadArchRProject(paste0(wd, "/filtered"))

##########################################################################################
# Copy ArchR project to new directory for adding peak information
##########################################################################################

proj <- saveArchRProject(
  ArchRProj = orig_proj,
  outputDirectory = outdir
)

proj <- loadArchRProject(outdir)

##########################################################################################
# Call Peaks
##########################################################################################

# Create Group Coverage Files that can be used for downstream analysis
proj <- addGroupCoverages(
  ArchRProj=proj, 
  groupBy= "FineNamedClust",
  #minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40, BOR = 50)
  )

# Find Path to Macs2 binary
pathToMacs2 <- findMacs2()

# Call Reproducible Peaks w/ Macs2
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "FineNamedClust", 
    peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
    pathToMacs2 = pathToMacs2,
    force = TRUE
)

# Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)

# Save project
saveArchRProject(proj)

##########################################################################################
# Motif annnotations and chromVAR deviation
##########################################################################################

# Load and convert PFMS from FigR paper into PWMs
cisbp_pwms <- readRDS("/oak/stanford/groups/wjg/skim/resources/MotifPWMs/FigR_cisbp_pwms/cisBP_human_pfms_2021.rds")
cisbp_pwms <- TFBSTools::toPWM(cisbp_pwms)

# Add motif annotations based on Buenrostro lab's cleaned up cisbp pwms
proj <- addMotifAnnotations(ArchRProj = proj, name = "Motif", motifPWMs = cisbp_pwms)

# Add background peaks
proj <- addBgdPeaks(proj)

#(WARNING: 2h+ use more cores and memory, if it fails go to 1.0.2 release of ArchR)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(proj)

##########################################################################################
# Identifying Marker Peaks
##########################################################################################
proj <- loadArchRProject(outdir)

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "FineNamedClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peak-Marker-Heatmap-FineNamedClust", width=10, height=15, ArchRProj=proj, addDOC=FALSE)

saveRDS(markerPeaks, paste0(proj@projectMetadata$outputDirectory, "/Plots/lda_markerPeaks_FineNamedClust.rds"))

##########################################################################################
# Add UMAP from Seurat
##########################################################################################

#proj <- loadArchRProject(outdir)
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# Grab umap coordinates and adjust colnames and cell barcodes to add to archr project
umap <- obj@reductions$umap@cell.embeddings
CBs <- rownames(umap)
str_sub(CBs, -19, -19) <- "#"
rownames(umap) <- CBs
umap <- data.frame(umap)
umap <- umap[which(rownames(umap) %in% getCellNames(proj)),]
colnames(umap) <- c("custom#UMAP1", "custom#UMAP2")

# Add seurat umap coordinates into ArchR
proj@embeddings$customUMAP <- SimpleList(df = umap, params = list())

# Save project
saveArchRProject(proj)

##########################################################################################
# Plot Motif Enrichment in Marker Peaks
##########################################################################################

#Identify Motif Enrichments in marker peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
    enrichMotifs, 
    n=5, 
    transpose=FALSE, 
    cutOff=log10pCut
)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-FineNamedClust", width=8, height=12, ArchRProj=proj, addDOC=FALSE)

saveArchRProject(proj)

##########################################################################################
# Plot motif deviations across all the clusters
##########################################################################################

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

##########################################################################################
# Add Seurat gene expression matrix to the ArchR Project 
##########################################################################################
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

proj <- loadArchRProject(outdir)

proj <- addSeuratGeneExpressionMatrix(archr.proj = proj, 
                                      gtf.file = "/oak/stanford/groups/wjg/skim/resources/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz",
                                      seu.file = obj, assay = "RNA", matchCellNames = TRUE)

saveArchRProject(proj)







