#!/usr/bin/env Rscript

#####################################
# Cluster scRNA using iterative LSI
#####################################

suppressPackageStartupMessages({
	library(dplyr)
	library(tidyr)
	library(Seurat)
	library(ggrastr)
	library(Rmagic)
	library(future) # For parallelization
	library(data.table)
})

#### Parameters ####
subgroup <- "preprocessed"
pointSize <- 0.2
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 10
plan("multicore", workers=nThreads)

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory and Plot directory to Folder
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/rna_preprocess_output"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

plotDir <- paste0(wd,"/clustering_qc_final")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# color palettes
sample_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/sample_cmap.rds")
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/final/gest_age_cmap.rds")

sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(obj$Sample)]
gest_age_cmap <- gest_age_cmap[names(gest_age_cmap) %in% unique(obj$age)]

##########################################
# Read in previously created Seurat object
##########################################
message("Reading in data...")
obj <- readRDS(paste0(wd, sprintf('/%s.rds', subgroup)))

DefaultAssay(obj) <- "SCT"

allGenes <- rownames(obj)


####################################################
# QC plots of clustering from previous seurat object
####################################################
# Plot clustering results:
message("Plotting clustering results...")
plotClusterQC(obj, subgroup="preprocessed", plotDir=plotDir, pointSize=pointSize, sampleCmap=sample_cmap, gest_age_cmap=gest_age_cmap)
message("Done.")


##########################################
# Identify markers per cluster
##########################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
message("Finding marker genes using Seurat...")
obj.markers <- FindAllMarkers(
    obj, 
    only.pos=TRUE, min.pct=0.25, logfc.threshold=1.0 # Don't bother using too many cells for this step
    )

# Save markers
write.table(obj.markers, file=paste0(wd, sprintf("/marker_genes_%s.tsv", subgroup)), 
    sep="\t", quote=FALSE,row.names=FALSE)


##########################################
# UMAP of high level groups
##########################################
message("Plotting selected marker features on UMAP...")

plotDir <- paste0(wd,"/expression_plots_final")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# Set colormaps
qualcmap <- cmaps_BOR$stallion
quantcmap <- cmaps_BOR$sunrise

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Epithelial" = c("EPCAM", "CDH1"),
    "Airway" = c("SOX2", "SCGB3A1" ,"SCGB3A2", "TPPP3"),
    "Ciliated" = c("FOXJ1", "CAPS"),
    "Basal" = c("TP63", "KRT5", "KRT17"),
    "PNECs" = c("ASCL1", "NEUROD1", "GRP", "GHRL"),
    "Alveolar" = c("SFTPB", "SFTPC", "AGER", "PDPN", "ETV5", "HHIP"),
    "Mesothelial" = c("MSLN", "C3"),
    "Endothelial" = c("PECAM1", "CLDN5", "CDH5"),
    "Venous" = c("ACKR1", "CPE"),
    "Aerocytes" = c("EDNRB", "S100A3"),
    "gCaps" = c("CA4", "IL7R"),
    "Lymphatic" = c("CCL21", "PROX1"),
    "Contractile" = c("ACTA2", "CNN1", "TAGLN"),
    "Fibroblasts" = c("COL1A1", "COL1A2", "PDGFRA", "DACH2"),
    "Chondrocytes" = c("ACAN", "COL2A1"),
    "AdvFibro" = c("SERPINF1", "COL3A1"),
    "AlvFibro" = c("MEOX2", "WNT2", "PLEKHH2"),
    "MyoF" = c("DACH2", "POSTN", "ROBO2", "EYA4"),
    "Peri" = c("TRPC6", "EBF1", "LRRTM4", "SULT1E1", "RBFOX1", "PAG1"),
    "VascSM" = c("SLIT3", "ELN", "ITGA11", "GPC3", "CSMD1", "PLN"),
    "AiwaySM" = c("HHIP", "HPSE2"),
    "NK" = c("NKG7", "GNLY"),
    "T_cells" = c("LTB", "CAMK4"),
    "B_cells" = c("IGHM", "TCF4"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "APCs" = c("LYZ", "CD86", "CD74") # Monocyte lineage (FCGR3A = CD16, FCGR1A = CD64, CD74 = HLA-DR antigens-associated invariant chain)
)

# Get expression data:
expr <- GetAssayData(obj, assay = "SCT", slot = "data") %>% t()
expr <- as(expr[,Matrix::colSums(expr) > 0], "sparseMatrix") # Remove unexpressed genes

selectedGenes <- unlist(featureSets) %>% unname()

flag <- "noMagic"
# Smooth with Rmagic
if(useMagic){
    message("Using MAGIC to impute (smooth) data for plotting...")

    # Run MAGIC directly on the expression matrix
    expr <- magic(expr, genes=selectedGenes, n.jobs = 1, seed = 1)$result
    flag <- "yesMagic"
}

for(name in names(featureSets)){
    features <- featureSets[[name]]
    pdf(paste0(plotDir,"/", name, "_features_UMAP.pdf"))
    for(gene in features){
        if(!gene %in% allGenes){
            message(sprintf("Error: %s is not a valid gene name", gene))
        }else{
            umapDF <- data.frame(Embeddings(object = obj, reduction = "umap"), expr[,gene])        
            colnames(umapDF) <- c("UMAP1", "UMAP2", gene)
            # Clip range of expression:
            upperLim <- quantile(umapDF$gene, probs=c(0.95))
            umapDF[,gene][umapDF[,gene] >= upperLim] <- upperLim
            print(plotUMAP(umapDF, dataType = "quantitative", cmap = quantcmap, covarLabel = gene, point_size=pointSize))
        } 
    }
    dev.off()
}

# Dot plot of cluster markers
count_mat <- GetAssayData(object = obj, assay = "SCT", slot = "counts")
avgPctMat <- avgAndPctExpressed(count_mat, obj$Clusters, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")

wide_df <- prettyOrderMat(wide_df)
grp_order <- colnames(wide_df$mat)
gene_order <- rev(rownames(wide_df$mat))

pdf(paste0(plotDir, "/markers_dot_plot_preclustered_final.pdf"), width=6, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise)
dev.off()




