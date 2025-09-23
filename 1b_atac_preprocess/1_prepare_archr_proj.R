#!/usr/bin/env Rscript

#####################################################################
# Build ArchR project and perform basic pre-processing and subsetting
#####################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(Seurat)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/sample_metadata.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess"

# color palettes
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds")
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations

pointSize <- 0.5
addArchRGenome("hg38")

##########################################################################################
# Preparing Data
##########################################################################################

#Get lung input files
inputFiles <- list.files("/oak/stanford/groups/wjg/skim/projects/LDA/0_cellranger/fragment_files", 
    recursive = TRUE, 
    full.names = TRUE,
    pattern = "*gz$" #ignores tsv.gz.tbi using the $
    ) 

names(inputFiles) <-  list.dirs("/oak/stanford/groups/wjg/skim/projects/LDA/0_cellranger/fragment_files/", 
    full.names = FALSE, 
    recursive = FALSE
    )

# Create Arrow Files (~30 minutes)
# recommend you use as many threads as samples.
# This step will for each sample :
# 1. Read Accessible Fragments.
# 2. Identify Cells QC Information (TSS Enrichment, Nucleosome info).
# 3. Filter Cells based on QC parameters.
# 4. Create a TileMatrix 500-bp genome-wide.
# 5. Create a GeneScoreMatrix.

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 5,
  minFrags = 1000, # Default is 1000.
  addTileMat = FALSE, # Don't add tile or geneScore matrices yet. Will add them after we filter
  addGeneScoreMat = FALSE
)

# If path to arrowfiles need to be redefined
#ArrowFiles <- list.files("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/unfiltered_output/ArrowFiles",
#                         full.names = TRUE)

#Create ArchRProject and save
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "raw_output"
)

saveArchRProject(proj)

# Remove Arrow files after copying
unlink(paste0(wd, "/*.arrow"))

proj <- loadArchRProject(paste0(wd, "/raw_output"))

# Filter cells based on final cell calls
# Grab Seurat Cell Barcodes and convert to ArchR compatible cell barcodes
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")
realCells <- colnames(obj)
stringr::str_sub(realCells, start = -19, end = -19) <- "#"
realCells <- realCells[which(realCells %in% getCellNames(proj))]

# Subset to real cells from RNA and ATAC doublet filtered cells
subProj <- subsetArchRProject(proj, cells=realCells, 
    outputDirectory="filtered", dropCells=TRUE)

# Add sample metadata
subProj$preservation <- samp.preservation[subProj$Sample] %>% unlist() %>% as.factor()
subProj$sex <- samp.sex[subProj$Sample] %>% unlist() %>% as.factor()
subProj$age <- samp.age[subProj$Sample] %>% unlist()

saveArchRProject(subProj, dropCells = TRUE)
subProj <- loadArchRProject(paste0(wd, "/filtered"))

# Now, add tile matrix and gene score matrix to ArchR project
subProj <- addTileMatrix(subProj, force=TRUE)
subProj <- addGeneScoreMatrix(subProj, force=TRUE)

# Visualize numeric metadata per grouping with a violin plot now that we have created an ArchR Project.
plotList <- list()

plotList[[1]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  pal=sample_cmap,
  plotAs = "violin",
  addBoxPlot = TRUE,
  alpha = 0.4
  )

plotList[[2]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  pal=sample_cmap,
  plotAs = "violin",
  addBoxPlot = TRUE,
  alpha = 0.4
  )

plotPDF(plotList = plotList, name = "Violin-TSS-log10(nFrag)", width = 4, height = 4,  ArchRProj = subProj, addDOC = FALSE)

# Save filtered ArchR project
saveArchRProject(subProj)

##########################################################################################
# Add gene expression matrix to ArchRProject
##########################################################################################

subProj <- addSeuratGeneExpressionMatrix(archr.proj = subProj,
                                         gtf.file = "/oak/stanford/groups/wjg/skim/resources/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz",
                                         seu.file = obj, matchCellNames = TRUE, assay = "RNA")

saveArchRProject(subProj)

##########################################################################################
# Reduced Dimensions and Clustering
##########################################################################################

# Reduce Dimensions with Iterative LSI (<5 minutes)
set.seed(1)
subProj <- addIterativeLSI(
    ArchRProj = subProj,
    useMatrix = "TileMatrix", 
    name = "LSI_ATAC", 
    sampleCellsPre = 25000,
    varFeatures = 50000, 
    dimsToUse = 1:30,
    force = TRUE
)

# Identify Clusters from Iterative LSI
subProj <- addClusters(
    input = subProj,
    reducedDims = "LSI_ATAC",
    method = "Seurat",
    name = "Clusters_ATAC",
    resolution = 0.6,
    force = TRUE
)

saveArchRProject(subProj)

##########################################################################################
# Incorporate RNA clustering results
##########################################################################################

subProj <- loadArchRProject(paste0(wd, "/filtered"))

# Add RNA clustering information from seurat object
rnaMetaData <- obj@meta.data
rnaMetaData$CB_Seurat <- rownames(rnaMetaData)
rnaMetaData$CB_ArchR <- rownames(rnaMetaData)
stringr::str_sub(rnaMetaData$CB_ArchR, start = -19, end = -19) <- "#"

atacMetaData <- data.frame(CB_ArchR = subProj$cellNames)
atacMetaData <- left_join(atacMetaData, rnaMetaData, by = "CB_ArchR")

# Add BroadNamed clusters from RNA
subProj <- addCellColData(
    ArchRProj = subProj,
    name = "BroadNamedClust",
    cells = atacMetaData$CB_ArchR,
    data = paste0(atacMetaData$BroadNamedClust),
    force = TRUE
    )

# Add FineNamed clusters from RNA
subProj <- addCellColData(
    ArchRProj = subProj,
    name = "FineNamedClust",
    cells = atacMetaData$CB_ArchR,
    data = paste0(atacMetaData$FineNamedClust),
    force = TRUE
    )

# Add fine cluster coded names from RNA
subProj <- addCellColData(
    ArchRProj = subProj,
    name = "FineClust",
    cells = atacMetaData$CB_ArchR,
    data = paste0(atacMetaData$FineClust),
    force = TRUE
    )

#Add Broad clusters from RNA
subProj <- addCellColData(
    ArchRProj = subProj,
    name = "BroadClust",
    cells = atacMetaData$CB_ArchR,
    data = paste0(atacMetaData$BroadClust),
    force = TRUE
    )

# Add Compartment information
subProj <- addCellColData(
    ArchRProj = subProj,
    name = "compartment",
    cells = atacMetaData$CB_ArchR,
    data = paste0(atacMetaData$compartment),
    force = TRUE
    )


##########################################################################################
# Visualize Data
##########################################################################################

sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(subProj$Sample)]
samp_cmap <- unlist(sample_cmap)

gest_age_cmap <- gest_age_cmap[names(gest_age_cmap) %in% unique(subProj$age)]
age_cmap <- unlist(gest_age_cmap)

set.seed(1)
subProj <- addUMAP(
    ArchRProj = subProj, 
    reducedDims = "LSI_ATAC", 
    name = "UMAP_ATAC", 
    nNeighbors = 50, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)

# Relabel clusters so they are sorted by cluster size
subProj <- relabelClusters(subProj, clusterName="Clusters_ATAC")

#subProj <- addImputeWeights(subProj)

# Make various cluster plots:
subProj <- visualizeClustering(subProj, 
    clusterName="BroadNamedClust", sampleName="Sample", embedding="UMAP_ATAC",
    pointSize=pointSize, sampleCmap=samp_cmap, gest_age_cmap=age_cmap,
    filename = "Plot-ATAC_UMAP-Sample-RNA_BroadNamedClusters.pdf")

subProj <- visualizeClustering(subProj, 
    clusterName="FineNamedClust", sampleName="Sample", embedding="UMAP_ATAC",
    pointSize=pointSize, sampleCmap=samp_cmap, gest_age_cmap=age_cmap,
    filename = "Plot-ATAC_UMAP-Sample-RNA_FineNamedClusters.pdf")

subProj <- visualizeClustering(subProj, 
    clusterName="Clusters_ATAC", sampleName="Sample", embedding="UMAP_ATAC",
    pointSize=pointSize, sampleCmap=samp_cmap, gest_age_cmap=age_cmap,
    filename = "Plot-ATAC_UMAP-Sample-ATAC_Clusters.pdf")

# Save filtered ArchR project
saveArchRProject(subProj)


##########################################################################################
# Plot cluster concordance between ATAC and RNA
##########################################################################################
proj <- subProj
proj <- loadArchRProject(paste0(wd, "/filtered_output/"), force=TRUE)

cM <- as.matrix(confusionMatrix(proj$BroadNamedClust, proj$Clusters_ATAC))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}
cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

pdf(paste0(wd, "/filtered_output/Plots/", "/ATAC-RNA-Multiome-cM-heatmap_BroadNamedClust.pdf"), width=6, height=6)
hm <- BORHeatmap(
  cM, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  row_names_side = "left",
  width = ncol(cM)*unit(0.5, "cm"),
  height = nrow(cM)*unit(0.5, "cm"),
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()
