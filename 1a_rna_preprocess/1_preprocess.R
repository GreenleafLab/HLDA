#!/share/software/user/open/R/4.0.2/bin/Rscript

#####################################
# Preprocess scRNA using Seurat
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
  library(Seurat)
  library(Matrix)
  library(stringr)
})

#### Parameters ####

# Cell quality filter criteria
minFeatures <- 200
maxFeatures <- Inf
minCounts <- 1000
maxCounts <- Inf
maxPctMito <- 20

# Seurat parameters
seuratFeatures <- 3000
seuratDims <- 1:30
seuratRes <- 0.4 # Lower resolution will increase estimated RNA contamination

# UMAP parameters
umapDims <- 40
umapNeighbors <- 50
umapMinDist <- 0.5
umapDistMetric <- "cosine"

# Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/rna_preprocess"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

# Set/Create plot directory
plotDir <- paste0(wd,"/raw_qc")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# Cellranger output data directory 
raw_data_dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/0_cellranger/filtered_feature_bc_matrix" # User defined


#### Logging ####
# Start logging
logfile <- paste0(wd, sprintf("/preprocess_log_%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")))
con <- file(logfile, open = "wt")
sink(con, type="output")
sink(con, type="message")

# Print all parameters to log file
for ( obj in ls() ) { cat('---',obj,'---\n'); print(get(obj)) }

#### Directories ####

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))

# Read in sample metadata
source(paste0(scriptPath, "/sample_metadata.R")) # Confirm metadata


#### RNA preprocess ####

# Load each of the datasets
message("Reading in RNA data...")
data_dirs <- list.dirs(path=raw_data_dir, full.names=TRUE, recursive=FALSE)
names(data_dirs) <- str_extract(data_dirs, "(?<=/)[^/]*$")

# Create individual Seurat Objects
objs <- list()
for(ix in seq_along(data_dirs)){
  sample <- names(data_dirs)[ix]
  path <- data_dirs[ix]
  message(sprintf("Reading in data from sample %s...", sample))
  data <- Read10X(data.dir=path)
  obj <- CreateSeuratObject(counts=data$`Gene Expression`, project=sample, min.cells=0, min.features=minFeatures)
  obj$orig.ident <- sample
  objs[[sample]] <- obj
}

# Extract names for each object
objNames <- sapply(objs, function(x) x@project.name)

# Close connections since individual seuratobject logging takes over
sink(type = "message")
sink(type = "output")
close(con)

# Perform initial scRNA preprocessing (cell calling, quality filtering, doublet detection, ambient RNA removal):
objs <- lapply(seq_along(objs), function(i){
  objName <- objNames[i]
  obj <- objs[[i]]
  scRNAdataPreProcessing(
    obj, objName, plotDir, # Naming parameters
    minFeatures=minFeatures, maxFeatures=maxFeatures, 
    minCounts=minCounts, maxCounts=maxCounts, maxPctMito=maxPctMito, # Quality filters
    nfeatures=seuratFeatures, dims=seuratDims, res=seuratRes, # Seurat clustering parameters
    runDoubletFinder=useDoubletFinder, runDecontX=useDecontX, # Optional processing steps
    estDubRate=estDubRate, # DoubletDetector parameters
    ncores=1
    )
  })

# Reopen main log
con <- file(logfile, open = "at")
sink(con, type="output", append=TRUE)
sink(con, type="message", append=TRUE)

# Save filtered Seurat object list:
message("Saving list of filtered individual Seurat object...")
saveRDS(objs, file = paste0(wd, "/preprocessed_list.rds"))
#objs <- readRDS(paste0(wd, "/preprocessed_list.rds"))

# Merge individual Seurat objects
obj <- merge(x=objs[[1]], y=objs[2:length(objs)], project="lung")

# Read in cell IDs passing filter
CellsPassingFilterDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/8-LDA/results/CellsPassingFilter.rds"
message(sprintf("Reading in final cell barcodes passing filter from %s", CellsPassingFilterDir))
CellsPassingFilter <- readRDS(CellsPassingFilterDir)

# Subset to final cells after manual doublet removal
cellsBeforeFiltering <- dim(obj)[2]
obj <- obj[, which(colnames(obj) %in% CellsPassingFilter)]
cellsAfterFiltering <- dim(obj)[2]
message(sprintf("%s filtered down to %s (%s%% remaining)", 
  cellsBeforeFiltering, cellsAfterFiltering, 
  round(100*(cellsAfterFiltering/cellsBeforeFiltering), 2)))

# Store sample and other metadata in object
obj$Sample <- obj$orig.ident
obj$preservation <- samp.preservation[obj$Sample] %>% unlist() %>% as.factor()
obj$sex <- samp.sex[obj$Sample] %>% unlist() %>% as.factor()
obj$age <- samp.age[obj$Sample] %>% unlist()

saveRDS(obj, file = paste0(wd, "/preprocessed_list_merged.rds"))

# Identify genes we want to blacklist from variable gene selection
# mitochondrial:
mt.genes <- grep(pattern = "^MT-", x = rownames(obj), value = TRUE)

# Cell cycle (These are loaded by Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# X/Y chromosome genes:
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

geneGR <- GenomicFeatures::genes(txdb)
sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
require(org.Hs.eg.db)
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
    g2m.genes
)


saveRDS(blacklist.genes, paste0(wd, "/blacklist.genes.rds"))

# Filtering the object list to remove doublets and low quality cells prior to variable feature calling
SampleNames <- obj$orig.ident %>% unique()

# Filter the list of objs based on cells that pass the final filter
objs_filtered <- lapply(seq_along(objs), function(i){
  SampleName <- SampleNames[i]
  so <- objs[[i]]
  so <- so[, which(colnames(so) %in% CellsPassingFilter)]
})

# Get consensus variable features across individually SCTransformed objects
# Based on: https://github.com/satijalab/seurat/issues/5761 with filtering from BOR
consFeatures <- getConsensusVarFeatures(objs_filtered, nfeatures = seuratFeatures, blacklist.genes = blacklist.genes)
saveRDS(consFeatures, paste0(wd, "/consensusFeatures.rds"))
#consFeatures <- readRDS(paste0(wd, "consensusFeatures.rds"))

# Fully scale data on merged object
obj <- SCTransform(obj, residual.features = consFeatures, vst.flavor = "v2")

# PCA, UMAP, Clustering
VariableFeatures(obj) <- consFeatures
obj <- RunPCA(obj, features = consFeatures) 
obj <- RunUMAP(obj, dims = 1:umapDims, n.neighbors = umapNeighbors, min.dist = umapMinDist) #umapDims of 40 will replicate the original UMAP
obj <- FindNeighbors(obj, dims = 1:umapDims)
obj <- FindClusters(obj)

pdf(paste0(wd, "/timepoint_UMAP_split_fullSCTransform.pdf"))
print(DimPlot(obj))
timepoints <- SplitObject(obj, split.by = "age")
for(i in timepoints) {print(DimPlot(i))}
dev.off()

# Store cluster information
obj$Clusters <- Idents(obj)

# Add cell cycle information now (Using Seurat):
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

# Depth normalize to 10,000, add pseudo count, and then log2 transform and add variable features
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
VariableFeatures(obj) <- consFeatures

saveRDS(obj, file = paste0(wd, "/preprocessed.rds"))


# Close connections
on.exit({ sink(type = "message"); sink(type = "output"); close(con) })

