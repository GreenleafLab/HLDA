#!/usr/bin/env Rscript

#####################################
# Cluster scRNA using iterative LSI
#####################################

# Subcluster previously identified groups

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(future)
  library(Matrix)
  library(harmony)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

#### Parameters ####

# Set/Create Working Directory to Folder:
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/rna_preprocess_output/subclustering_final"
setwd(wd)

# Subgroups to subcluster:
subgroups <- c("Epithelial", "Endothelial", "Stromal", "Immune")

# Subclustering parameters:
paramDict <- list(
  "Immune" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "clusterRes" = 0.2,
    "nNeighbors" = 50,
    "minDist" = 0.5,
    "pointSize" = 0.25
    ),
  "Epithelial" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "clusterRes" = 0.2,
    "nNeighbors" = 50,
    "minDist" = 0.5,
    "pointSize" = 0.25
    ),
  "Stromal" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "clusterRes" = 0.2,
    "nNeighbors" = 50,
    "minDist" = 0.5,
    "pointSize" = 0.25
    ),
  "Endothelial" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "clusterRes" = 0.2,
    "nNeighbors" = 50,
    "minDist" = 0.5,
    "pointSize" = 0.25
    )
  )

# Misc
#nThreads <- 10
umapDistMetric <- "cosine"

# change the current plan to access parallelization (for Seurat)
#plan("multicore", workers = nThreads)

# Start logging:
logfile <- paste0(wd, sprintf("/subclustering_log_%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")))
con <- file(logfile, open = "wt")
sink(con, type="output")
sink(con, type="message")

# Print all parameters to log file
for ( obj in ls() ) { cat('---',obj,'---\n'); print(get(obj)) }

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/iterative_LSI.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# color palettes
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/sample_cmap.rds")
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/gest_age_cmap.rds")


# Identify genes we want to blacklist during clustering

# First get list of all genes:
subwd <- paste0(wd, sprintf("/%s", subgroups[1]))
allGenes <- rownames(GetAssayData(object=readRDS(paste0(subwd, sprintf('/%s.rds', subgroups[1]))), slot="counts"))

# Blacklist genes
# mitochondrial:
mt.genes <- grep(pattern="^MT-", x=allGenes, value=TRUE)
# Cell cycle: (These are loaded by Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# Ribosomal:
rp.genes <- grep(pattern="^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", x=allGenes, value=TRUE)

# X/Y chromosome genes:
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
geneGR <- GenomicFeatures::genes(txdb)
sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
matchedGeneSymbols <- select(org.Hs.eg.db,
                        keys = sexGenesGR$gene_id,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "ENTREZID")
sexChr.genes <- matchedGeneSymbols$SYMBOL


# Genes to ignore (just for clustering purposes)
blacklist.genes <- c(
    mt.genes,
    sexChr.genes,
    s.genes,
    g2m.genes,
    rp.genes
)



# Iteratively cluster subgroups:
for(sg in subgroups){

  subwd <- paste0(wd, sprintf("/%s", sg))

  # Directory for clustering qc plots:
  plotDir <- paste0(subwd, "/clustering_qc")
  dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

  # Read in previously created Seurat subobjects:
  message(sprintf("Reading in data for subgroup %s...", sg))
  obj <- readRDS(paste0(subwd, sprintf('/%s.rds', sg)))

  # # Remove any existing DimReductions
  # obj <- DietSeurat(obj, features=NULL, assays=NULL, dimreducs=NULL)

  # Subset colors to only those samples present
  samp_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)] %>% unlist()

  #######################################
  # Perform clustering with consensus feature set
  #######################################

  # Get subgroup parameters:
  nVarGenes <- paramDict[[sg]]$nVarGenes
  nPCs <- paramDict[[sg]]$nPCs
  clusterRes <- paramDict[[sg]]$clusterRes

  # UMAP:
  umapNeighbors <- paramDict[[sg]]$nNeighbors
  umapMinDist <- paramDict[[sg]]$minDist

  # Split obj based on sample
  objs <- SplitObject(obj, split.by = "Sample")

  # Perform sctransform for each sample
  objs <- lapply(X = objs, FUN = function(x) {
    x <- SCTransform(x, vst.flavor = "v2")
    })

  # Identify consensus features across all the samples
  message("Identifying consensus variable features")
  consFeatures <- getConsensusVarFeatures(objs, nfeatures = nVarGenes, blacklist.genes = blacklist.genes)

  # Merge seurat objects into one with renormalization using the residuals
  obj <- merge(x=objs[[1]], y=objs[2:length(objs)], project=sg)

  obj <- SCTransform(obj, residual.features = consFeatures, vst.flavor = "v2")

  # PCA, UMAP, Clustering
  obj <- RunPCA(obj, features = consFeatures) 

  message("Calculating UMAP...")
  obj <- RunUMAP(obj, dims = nPCs, n.neighbors = umapNeighbors, min.dist = umapMinDist)
  obj <- FindNeighbors(obj, dims = nPCs)
  obj <- FindClusters(obj, resolution = clusterRes)

  # Store cluster information in metadata
  obj$Clusters <- Idents(obj)

  message("Saving seurat object...")
  # Save clustered object here:
  saveRDS(obj, file = paste0(subwd, sprintf("/%s.rds", sg)))

  ##################################################
  # Plot Clustering results
  ##################################################
  # Plot clustering results:
  message("Plotting clustering results...")
  pointSize <- paramDict[[sg]]$pointSize
  plotClusterQC(obj, subgroup=sg, plotDir=plotDir, pointSize=pointSize, sampleCmap=sample_cmap, gest_age_cmap=gest_age_cmap)

  message("Done.")

}



