#!/usr/bin/env Rscript

############################################
# Identify marker genes and make some plots
############################################

library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr)
library(Rmagic)
library(future)
library(data.table)

subgroup <- "lda"
pointSize <- 0.5
useMagic <- TRUE # Should Rmagic be used for data imputation prior to UMAP plotting?

# change the current plan to access parallelization (for Seurat)
nThreads <- 10
plan("multicore", workers = nThreads)

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Setup working directory and make a plot dir

#Set/Create Working Directory to Folder
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/final/rna_preprocess_output"
plotDir <- paste0(wd,"/expression_plots_final")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

##########################################
# Read in previously created Seurat object
##########################################
message("Reading in data...")
obj <- readRDS(paste0(wd, sprintf('/%s.rds', subgroup)))

DefaultAssay(obj) <- "SCT"

Idents(obj) <- "FineNamedClust"

allGenes <- rownames(obj)

##########################################
# Identify markers per cluster (And GO terms)
##########################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
message("Finding marker genes using Seurat...")
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)

# Save markers
write.table(obj.markers, file=paste0(wd, sprintf("/marker_genes_%s.tsv", subgroup)), 
    sep="\t", quote=FALSE,row.names=FALSE)

# Use only expressed genes in GO term analyses:
# Define 'expressed genes' as those with at least 2 counts in at least 5 cells
rawCounts <- GetAssayData(object = obj, slot = "counts")
minUMIs <- 2
minCells <- 5
expressedGenes <- rownames(rawCounts[rowSums(rawCounts > minUMIs) > minCells,])

rm(rawCounts); gc()

message("Calculating GO terms on each cluster using marker genes...")

clusters <- unique(obj.markers$cluster)
Idents(obj) <- "Clusters"

# Get GO enrichments (expressed genes):
GOresults <- lapply(clusters, function(x){
  message(sprintf("Running GO enrichments on cluster %s...", x))
  df <- obj.markers[obj.markers$cluster == x,]
  # Define custom cutoff for 'interesting genes'
  intGenes <- df[(df$p_val_adj < 0.01) & (df$avg_log2FC > 0.5),"gene"]
  calcTopGo(expressedGenes, interestingGenes=intGenes)
  })
names(GOresults) <- paste0("cluster", clusters)

# Plots of GO term enrichments:
pdf(paste0(plotDir, "/cluster_marker_GO_term_enrichments.pdf"), width=8, height=8)
for(name in names(GOresults)){
    goRes <- GOresults[[name]]
    print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
        nterms=10, border_color="black", 
        barwidth=0.85, title=name))
}
dev.off()

# Save a file of the results for each cluster:
goDir <- paste0(wd,"/go_term_enrichments")
dir.create(goDir, showWarnings = FALSE, recursive = TRUE)
for(n in names(GOresults)){
    write.table(GOresults[[n]], file=paste0(goDir, "/", n, "_go_terms.tsv"), quote=FALSE, sep='\t')
}

##########################################
# UMAPs and DotPlots of Known Marker Genes
##########################################

# Set colormaps
qualcmap <- cmaps_BOR$stallion
quantcmap <- cmaps_BOR$sunrise

# Markers for identifying broad classes of cells:
featureSets <- list(
    "Epithelial" = c("EPCAM", "CDH1"),
    "Airway" = c("SOX2", "SCGB3A2", "TPPP3"),
    "Ciliated" = c("FOXJ1", "CAPS"),
    "Basal" = c("TP63", "KRT5", "KRT17"),
    "PNECs" = c("ASCL1", "NEUROD1", "GRP", "GHRL", "P3H2", "PCDH7", "CDH19"),
    "Alveolar" = c("SFTPB", "SFTPC", "AGER", "PDPN", "ETV5", "HHIP"),
    "Mesothelial" = c("MSLN", "C3"),
    "Endothelial" = c("PECAM1", "CLDN5", "CDH5"),
    "Venous" = c("ACKR1", "CPE"),
    "Arteries" = c("SERPINE2", "FBLN5", "GJA4", "DKK2"),
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
    "B_cells" = c("IGHM", "TCF4", "BANK1", "JCHAIN"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "APCs" = c("LYZ", "CD86", "CD74"), # Monocyte lineage (FCGR3A = CD16, FCGR1A = CD64, CD74 = HLA-DR antigens-associated invariant chain)
    "Monocytes" = c("VCAN", "AQP9"),
    "DendriticCell" = c("HLA-DPB1", "HLA-DPA1", "HLA-DRB1"),
    "Neutrophil" = c("LTF", "DEFA3"),
    "Cycling" = c("TOP2A", "CENPH", "MKI67")
)

# Get expression data:
expr <- GetAssayData(obj, slot = 'data') %>% t()
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
count_mat <- GetAssayData(object = obj, slot = "counts")
avgPctMat <- avgAndPctExpressed(count_mat, obj$FineNamedClust, feature_normalize=TRUE, min_pct=5)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")

#wide_df <- prettyOrderMat(wide_df[,rnaOrder], clusterCols=FALSE)
wide_df <- prettyOrderMat(wide_df, clusterCols=TRUE)

grp_order <- colnames(wide_df$mat)
gene_order <- rev(rownames(wide_df$mat))

pdf(paste0(plotDir, "/markers_dot_plot_lda.pdf"), width=6, height=10)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", xorder=grp_order, yorder=gene_order, cmap=cmaps_BOR$sunrise)
dev.off()

