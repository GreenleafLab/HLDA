#/share/software/user/open/R/4.0.2/bin/Rscript
#sbatch -p wjg,biochem,sfgf -t 03:00:00 --mem=50G --cpus-per-task=8 --job-name=3_assign_RNA_clusters --output=assign_RNA.out --error=assign_RNA.err --wrap "Rscript scripts/rna_preprocess/3_assign_RNA_clusters"

#####################################################
# Subcluster scRNA 
#####################################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)
library(future) # For parallelization

# change the current plan to access parallelization (for Seurat)
#nThreads <- 10
#plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Setup working directory and make a plot dir
# Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/rna_preprocess_output"
plotDir <- paste0(wd,"/expression_plots_final")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# color palettes
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/sample_cmap.rds")

##########################################
# Read in previously created Seurat object
##########################################
message("Reading in data...")
obj <- readRDS(paste0(wd, '/preprocessed_final.rds'))

# prep sample color map
sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)]

allGenes <- rownames(obj)

pointSize <- 0.15

##########################################
# Annotate broad clusters without any subclustering
##########################################

# Store the original clustering as BroadClust
obj$BroadClust <- obj$Clusters

BroadNamedClust <- list(
    "0" = "gCap", # gCap
    "1" = "AT1l", # AT1l_1
    "2" = "aSMC", #
    "3" = "TiP",  # TiP1
    "4" = "Peri", #
    "5" = "AlvF", # 
    "6" = "AT1l", #
    "7" = "APr", # APr4
    "8" = "TiP", # TiP2
    "9" = "AT1l", #
    "10" = "APr", # APr2
    "11" = "AlvF", #
    "12" = "APCs", # 
    "13" = "EpiC", #EpiC = cycling epithelial cells
    "14" = "NK", #
    "15" = "AT2l", #
    "16" = "vSMC", #CSMD1+, PLN+, RCAN2+ vSMC (vSMC2)
    "17" = "Aero", # Aero = Aerocytes
    "18" = "PNEC", # 
    "19" = "Cili", #
    "20" = "APr", # APr3
    "21" = "MyoF", #
    "22" = "Artr", #
    "23" = "Chdr", #
    "24" = "AdvF", #
    "25" = "Tc", #
    "26" = "Veno", #
    "27" = "APr", # APr1
    "28" = "Lymp", #
    "29" = "Cili", #
    "30" = "Meso" #
)

obj$BroadNamedClust <- BroadNamedClust[obj$BroadClust] %>% unlist() %>% unname

##########################################
# Additional subclustering on specific clusters
##########################################

# Adding subclustering for vSMC cluster (16)
obj <- FindSubCluster(obj, cluster = "16", graph.name = "SCT_snn", resolution = 0.2) #16_0 and 16_1

obj$Clusters <- obj$sub.cluster
Idents(obj) <- "sub.cluster"

# Adding subclustering for pericytes (4)
obj <- FindSubCluster(obj, cluster = "4", graph.name = "SCT_snn", resolution = 0.1) #4_0 and 4_1

obj$Clusters <- obj$sub.cluster
Idents(obj) <- "sub.cluster"

# Adding subclustering for PNEC (18)
obj <- FindSubCluster(obj, cluster = "18", graph.name = "SCT_snn", resolution = 0.2) #18_0, 18_1, 18_2, 18_3

obj$Clusters <- obj$sub.cluster
Idents(obj) <- "sub.cluster"

# Adding subclustering for gCaps (0)
obj <- FindSubCluster(obj, cluster = "0", graph.name = "SCT_snn", resolution = 0.13) #0_0, 0_1

obj$Clusters <- obj$sub.cluster
Idents(obj) <- "sub.cluster"

# Adding subclustering for Aero (17)
obj <- FindSubCluster(obj, cluster = "17", graph.name = "SCT_snn", resolution = 0.1) #17_0, 17_1

obj$Clusters <- obj$sub.cluster
Idents(obj) <- "sub.cluster"

# Adding subclustering for Arteries (22)
obj <- FindSubCluster(obj, cluster = "22", graph.name = "SCT_snn", resolution = 0.2) #22_0, 22_1

obj$Clusters <- obj$sub.cluster
Idents(obj) <- "sub.cluster"

# Adding subclustering for APCs (12)
obj <- FindSubCluster(obj, cluster = "12", graph.name = "SCT_snn", resolution = 0.15) 
#12_0 (Monocyte VCAN+, AQP9+), 
#12_1 (Dendritic cell HLA-DPB1, DPA1 DRB1 HLA positive)
#12_2 (Bcells JCHAIN+, IGHM+ BANK1+ )
#12_3 (Neutrophil LTF+, DEFA3+ (secondary granules of neutrophils and cytotoxic peptides in neutrophil granules) )

obj$Clusters <- obj$sub.cluster

Idents(obj) <- "Clusters"

##########################################
# Assigning cell type identity to clusters
##########################################
FineNamedClust <- list(
    "0_0" = "egCap", # early gCap
    "0_1" = "lgCap", # Late gCap
    "1" = "AT1l1", # 
    "2" = "aSMC", #
    "3" = "TiP1", #
    "4_0" = "ePeri", #SULT1E1+ pericytes (ePeri)
    "4_1" = "lPeri", #LRRTM4+, RBFOX1+, PAG1+, pericytes (lPeri)
    "5" = "lAlvF", # 
    "6" = "AT1l2", #
    "7" = "APr4", #
    "8" = "TiP2", #
    "9" = "AT1l3", #
    "10" = "APr2", #
    "11" = "eAlvF", #
    "12_0" = "Mono", # (Monocyte VCAN+, AQP9+), 
    "12_1" = "Dc", #(Dendritic cell HLA-DPB1, DPA1 DRB1 HLA positive)
    "12_2" = "Bc", #(Bcells JCHAIN+, IGHM+ BANK1+ )
    "12_3" = "Neut", # Neutrophils LTF+, DEFA3+ (secondary granules of neutrophils and cytotoxic peptides in neutrophil granules) 
    "13" = "EpiC", #EpiC = cycling epithelial cells
    "14" = "NK", #
    "15" = "AT2l", #
    "16_0" = "vSMC1", 
    "16_1" = "vSMC2", #CSMD1+, PLN+, RCAN2+ vSMC (vSMC2)
    "17_0" = "lAero", # Aero = Aerocytes
    "17_1" = "eAero",
    "18_0" = "PNEC2", # GRP+ PNEC
    "18_1" = "PNEC3", # GHRL+ PNEC
    "18_2" = "PNEC1", # ASCL1+ NEUROD1+
    "18_3" = "Schw", # CDH19+ Schwann
    "19" = "lCili", #
    "20" = "APr3", #
    "21" = "MyoF", #
    "22_0" = "eArtr", #
    "22_1" = "lArtr", #arteries SERPINE2+, FBLN5+
    "23" = "Chdr", #
    "24" = "AdvF", #
    "25" = "Tc", #
    "26" = "Veno", #
    "27" = "APr1", #
    "28" = "Lymp", #
    "29" = "eCili", #
    "30" = "Meso" #
)

# Store the fine cluster numbering
obj$FineClust <- obj$Clusters

# Add fine cluster names
obj$FineNamedClust <- FineNamedClust[obj$FineClust] %>% unlist() %>% unname

# Assign compartment information
compartments <- list(
  "Epithelial" = c("APr", "PNEC", "Cili", "TiP",  "AT2l", "AT1l", "EpiC", "Meso"),
  "Endothelial" = c("gCap", "Aero", "Artr", "Veno", "Lymp"),
  "Stromal" = c("AlvF", "AdvF", "MyoF", "aSMC", "Peri", "vSMC", "Chdr"),
  "Immune" = c("Tc", "NK", "APCs")
  ) %>% invertList()

# Unlist for adding info to seurat object
obj$compartment <- compartments[obj$BroadNamedClust] %>% unlist()

# Storing label hierarchy for color mapping
BroadlabelHierarchy <- list(
  "Epithelial" = c("APr", "PNEC", "Cili", "TiP",  "AT2l", "AT1l", "EpiC", "Meso"),
  "Endothelial" = c("gCap", "Aero", "Artr", "Veno", "Lymp"),
  "Stromal" = c("AlvF", "AdvF", "MyoF", "aSMC", "Peri", "vSMC", "Chdr"),
  "Immune" = c("Tc", "NK", "APCs")
  ) 

FinelabelHierarchy <- list(
  "Epithelial" = c("APr1", "APr2", "APr3", "APr4", "PNEC1", "PNEC2", "PNEC3", "eCili", "lCili", 
    "TiP1", "TiP2", "AT2l", "AT1l", "EpiC",  "Schw", "Meso"),
  "Endothelial" = c("egCap", "lgCap", "eAero", "lAero", "eArtr", "lArtr", "Veno", "Lymp"),
  "Stromal" = c("eAlvF", "lAlvF", "AdvF", "MyoF", "aSMC", "ePeri", "lPeri", "vSMC1", "vSMC2", "Chdr"),
  "Immune" = c("IM", "Mono1", "Mono2", "MPP", "Bc", "Neut", "Dc", "NK", "Plasma", "Tc")
  ) 

# Create color palette that will be used for all clustering
lungClusterColors <- c(
  "Immune" = "#D51F26", # red
  "Stromal" = "#272E6A", # dark blue
  "Epithelial" = "#208A42", # green
  "Endothelial" = "#F47D2B" # orange
  )

# Save these colors for use in other scripts
saveRDS(lungClusterColors, file = "/oak/stanford/groups/wjg/skim/projects/LDA/final/lungClusterColors.rds")
colorPal <- lungClusterColors
expandedColors <- getColorMap(cmaps_BOR$stallion, n=50)
expandedColors <- expandedColors[expandedColors %ni% colorPal]

# Now identify colormap for broadnamedclusters
broadLabels <- unique(obj$compartment)
narrowLabels <- unique(obj$BroadNamedClust)

narrowColors <- list()
takenColors <- colorPal
for(bl in broadLabels){
    subLabs <- BroadlabelHierarchy[[bl]]
    if(length(subLabs) == 1){
        # If only a single label, use broad label
        narrowColors[[subLabs[1]]] <- colorPal[[bl]]
    }else{
        validColors <- expandedColors[expandedColors %ni% takenColors]
        subCols <- mostSimilarColors(colorPal[[bl]], colorOptions=validColors, n=length(subLabs))
        for(i in seq_along(subLabs)){
            narrowColors[[subLabs[i]]] <- subCols[i]
            takenColors <- c(takenColors, subCols[i])
        }
    }
}
# Save color palette for 'BroadNamedClust'
saveRDS(narrowColors, file = "/oak/stanford/groups/wjg/skim/projects/LDA/final/scRNA_BroadNamedClust_cmap.rds")

# Now identify colormap for finenamedclusters
broadLabels <- unique(obj$compartment)
narrowLabels <- unique(obj$FineNamedClust)

narrowColors <- list()
takenColors <- colorPal
for(bl in broadLabels){
    subLabs <- FinelabelHierarchy[[bl]]
    if(length(subLabs) == 1){
        # If only a single label, use broad label
        narrowColors[[subLabs[1]]] <- colorPal[[bl]]
    }else{
        validColors <- expandedColors[expandedColors %ni% takenColors]
        subCols <- mostSimilarColors(colorPal[[bl]], colorOptions=validColors, n=length(subLabs))
        for(i in seq_along(subLabs)){
            narrowColors[[subLabs[i]]] <- subCols[i]
            takenColors <- c(takenColors, subCols[i])
        }
    }
}
# Save color palette for 'FineNamedClust'
saveRDS(narrowColors, file = "/oak/stanford/groups/wjg/skim/projects/LDA/final/scRNA_FineNamedClust_cmap.rds")

# Load colormaps for plotting:
barwidth=0.9
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/sample_cmap.rds")
sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)] %>% unlist()
lungClusterColors <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/lungClusterColors.rds") %>% unlist()
broadColors <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/scRNA_BroadNamedClust_cmap.rds") %>% unlist()
narrowColors <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/scRNA_FineNamedClust_cmap.rds") %>% unlist()

qualcmap <- cmaps_BOR$stallion

### Fine Named cluster UMAP ###
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), obj$FineNamedClust)
# Randomize cells before plotting
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir,"/FineNamedClusters_UMAP.pdf"), width=10, height=10)
plotUMAP(umapDF, dataType="qualitative", cmap=narrowColors, namedColors=TRUE, point_size=pointSize)
dev.off()

clustBySamp <- fractionXbyY(obj$FineNamedClust, obj$Sample, add_total=TRUE, xname="FineNamedClust", yname="Sample")
pdf(paste0(plotDir, "/FineclustBySampleBarPlot_NamedClusters.pdf"))
print(stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()


gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/gest_age_cmap.rds")
gest_age_cmap <- gest_age_cmap[names(gest_age_cmap) %in% unique(obj$age)]
clustByAge <- fractionXbyY(obj$FineNamedClust, obj$age, add_total=TRUE, xname="FineNamedClust", yname="Gestational Age")
pdf(paste0(plotDir, "/FineclustByAgeBarPlot.pdf"))
print(stackedBarPlot(clustByAge, cmap=gest_age_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

### Broad Named cluster UMAP ###
umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), obj$BroadNamedClust)
# Randomize cells before plotting
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir,"/BroadNamedClusters_UMAP.pdf"), width=10, height=10)
plotUMAP(umapDF, dataType="qualitative", cmap=broadColors, namedColors=TRUE, point_size=pointSize)
dev.off()

clustBySamp <- fractionXbyY(obj$BroadNamedClust, obj$Sample, add_total=TRUE, xname="BroadNamedClust", yname="Sample")
pdf(paste0(plotDir, "/BroadclustBySampleBarPlot_NamedClusters.pdf"))
print(stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()


gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/gest_age_cmap.rds")
gest_age_cmap <- gest_age_cmap[names(gest_age_cmap) %in% unique(obj$age)]
clustByAge <- fractionXbyY(obj$BroadNamedClust, obj$age, add_total=TRUE, xname="BroadNamedClust", yname="Gestational Age")
pdf(paste0(plotDir, "/BroadclustByAgeBarPlot.pdf"))
print(stackedBarPlot(clustByAge, cmap=gest_age_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

### Compartment to cluster figures ###
umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$compartment)
# Randomize cells before plotting
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]

pdf(paste0(plotDir,"/Compartments_UMAP.pdf"))
plotUMAP(umapDF, dataType="qualitative", namedColors=TRUE, point_size=pointSize)
dev.off()

clustBySamp <- fractionXbyY(obj$compartment, obj$Sample, add_total=TRUE, xname="Compartment", yname="Sample")
qualcmap <- cmaps_BOR$stallion
pdf(paste0(plotDir, "/clustBySampleBarPlot_Compartment.pdf"))
print(stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

clustByAge <- fractionXbyY(obj$compartment, obj$age, add_total=TRUE, xname="Compartment", yname="Gestational Age")
pdf(paste0(plotDir, "/clustByAgeBarPlot_BroadClusters.pdf"))
print(stackedBarPlot(clustByAge, cmap=gest_age_cmap, namedColors=TRUE, barwidth=barwidth))
dev.off()

# Save whole project with all cluster information:
saveRDS(obj, file = paste0(wd, "/lda.rds"))

########################################################################################################
# Subcluster groups of interest
########################################################################################################

obj <- readRDS(paste0(wd, "/lda.rds"))

# Make new Seurat objects for each sub-clustered group
DefaultAssay(obj) <- "RNA"

makeSubClusts <- function(obj, ident, subgroups, outdir){
  Idents(obj) <- ident
  for(subg in subgroups){
    subsubdir <- paste0(outdir, sprintf("/%s", subg))
    dir.create(subsubdir, showWarnings = FALSE, recursive = TRUE)
    subObj <- subset(obj, idents = c(subg))
    counts <- GetAssayData(object = subObj, slot = "counts")
    newObj <- CreateSeuratObject(counts = counts, project = subg, min.cells = 0, min.features = 200)
    old.meta <- subObj@meta.data
    # Drop selected columns from old metadata
    old.cols <- colnames(old.meta)
    drop.cols <- old.cols[grepl("^SCT_snn", old.cols)]
    newObj@meta.data <- old.meta[,old.cols %ni% drop.cols]
    message(sprintf("Subcluster %s has %s cells", subg, dim(newObj)[2]))
    saveRDS(newObj, file = paste0(subsubdir, "/", subg, ".rds"))
  }
}

subclustDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/rna_preprocess_output/subclustering_final"
dir.create(subclustDir, showWarnings = FALSE, recursive = TRUE)

makeSubClusts(
  obj, 
  ident="compartment", 
  subgroups=c("Epithelial", "Endothelial", "Stromal", "Immune"),
  outdir=subclustDir
)
