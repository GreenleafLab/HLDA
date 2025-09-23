#!/usr/bin/env Rscript

########################################
# Figure 2 panels
########################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(BuenColors)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/Figure02_Endothelial"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

# set plot directory
plotDir <- paste0(wd, "/plots")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# Misc options
addArchRGenome("hg38")
pointSize <- 1

##########################################################################################
# Load objects
##########################################################################################
# Load data objects
atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_v2", force = TRUE)
rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# Color Maps
compartmentCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/lungClusterColors.rds") %>% unlist()
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds") %>% unlist()
sample_cmap <- readRDS(paste0("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds"))
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")

# Load subprojects
sub_rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Endothelial/Endothelial.rds")
sub_atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Endothelial")

Idents(sub_rna_proj) <- "FineNamedClust"

##########################################################################################
# UMAPs
##########################################################################################
umapPlots <- list()

plot.name <- "UMAP_RNA_FineNamedClust"
dir.name <- paste0(wd, "/", plot.name)
dir.create(dir.name, showWarnings = FALSE, recursive = TRUE)

# FineNamedClust on RNA:
umapDF <- data.frame(Embeddings(object = sub_rna_proj, reduction = "umap"), sub_rna_proj$FineNamedClust)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
readr::write_tsv(umapDF, file = paste0(dir.name, "/UMAP_RNA_FineNamedClust.tsv"))
umapPlots[["RNA_FineNamedClust"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=FineNamedClustCmap, 
         namedColors=TRUE, point_size=pointSize, covarLabel="FineNamedClust_on_RNA", useRaster=TRUE)

# Age on RNA:
umapDF <- data.frame(Embeddings(object = sub_rna_proj, reduction = "umap"), sub_rna_proj$age)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
readr::write_tsv(umapDF, file = paste0(dir.name, "/UMAP_RNA_age.tsv"))
umapPlots[["RNA_age"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=gest_age_cmap, 
                                              namedColors=TRUE, point_size=pointSize, covarLabel="Gestational_age", useRaster=TRUE)

pdf(paste0(plotDir,"/",plot.name,".pdf"), width=7, height=5)
umapPlots
dev.off()

##########################################################################################
# Differential abundance testing for endothelial cell types
##########################################################################################
library(miloR)

#Set/Create Working Directory to Folder
subPlotDir <- paste0(wd, "/milo")
dir.create(subPlotDir, showWarnings = FALSE, recursive = TRUE)

Idents(rna_proj) <- "compartment"
endo_proj <- subset(rna_proj, ident = "Endothelial")

endo_proj@reductions$umap <- sub_rna_proj@reductions$umap

# Convert to SCE since milo needs it
sce <- as.SingleCellExperiment(endo_proj, assay = "SCT")
sce <- Milo(sce)

# Construct KNN graph
d <- 50
k <- 50

sce <- buildGraph(sce, k = k, d = d, reduced.dim = "PCA")

# Define representative neighborhoods
prop <- 0.1 # "We suggest using prop=0.1 for datasets of less than 30k cells. For bigger datasets using prop=0.05 should be sufficient (and makes computation faster)"
sce <- makeNhoods(sce, prop=prop, k=k, d=d, refined=TRUE, reduced_dims= "PCA")

# Plot the distribution of neighborhood sizes
pdf(paste0(subPlotDir, "/miloR_nhood_size_hist.pdf"), width = 8, height = 6)
plotNhoodSizeHist(sce)
#As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples
dev.off()

# Counting cells in neighborhoods
sce <- countCells(sce, meta.data = as.data.frame(colData(sce)), sample = "Sample")
head(nhoodCounts(sce))

# Define experimental design
design <- data.frame(colData(sce))[,c("Sample", "age")]
design <- distinct(design)
rownames(design) <- design$Sample

# Compute neighborhood connectivity (longest step)
sce <- calcNhoodDistance(sce, d=d, reduced.dim = "PCA")

# Testing differential abudance
da_results <- testNhoods(sce, design = ~ age, design.df = design, reduced.dim = "PCA")

# Visualize results:
sce <- buildNhoodGraph(sce)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce, da_results, layout="UMAP",
                                alpha=0.1, size_range=c(1, 4), node_stroke=0.1)

pdf(paste0(subPlotDir, "/miloR_RNA_DA_UMAP.pdf"), width=7, height=6)
nh_graph_pl
dev.off()

# Plot DA fold changes in different cell types
# Label neighborhoods with the cell type names
da_results <- annotateNhoods(sce, da_results, coldata_col = "FineNamedClust")
da_results$FineNamedClust <- ifelse(da_results$FineNamedClust_fraction < 0.7, "Mixed", da_results$FineNamedClust)
saveRDS(da_results, paste0(subPlotDir, "/milo_da_results.rds"))

# Check distribution of neighborhoods that are a mix of cell types
ggplot(da_results, aes(FineNamedClust_fraction)) + geom_histogram(bins = 50)

# Reorder the cell types
endo.order <- c(unlist(unname(fineOrder))[which(unlist(unname(fineOrder)) %in% unique(sub_rna_proj$FineNamedClust))], "Mixed")
da_results$FineNamedClust <- factor(da_results$FineNamedClust, levels = rev(endo.order))

# Plot beeswarm
pdf(paste0(subPlotDir, "/miloR_RNA_DA_beeswarm.pdf"), w = 5, h = 3)
plotDAbeeswarm(da_results, group.by = "FineNamedClust") + 
  guides(color=guide_colorbar()) +
  theme_BOR() +
  geom_hline(yintercept = 0)
dev.off()

##########################################################################################
# Fate probabilities
##########################################################################################

# Refer to Jupyter notebook 3_trajectories/Endothelial/01_WOT_Endothelial.ipynb

##########################################################################################
# Major trajectories
##########################################################################################
# Trajectory directory
trajDir <- paste0(wd, "/trajectories")
dir.create(trajDir, showWarnings = FALSE, recursive = TRUE)

# Endothelial trajectory (egCap -> eAero -> lAero)
# Define the endothelial trajectory based on optimal transport results
trajectory <- c("egCap", "eAero", "lAero")
trajectory_name <- "egCap_eAero_lAero"

sub_atac_proj <- addTrajectory(ArchRProj = sub_atac_proj, 
                               name = trajectory_name, 
                               groupBy = "FineNamedClust", 
                               trajectory = trajectory, 
                               reducedDims = "PCA",
                               force = TRUE
                               #embedding = "customUMAP"
)

p <- plotTrajectory(sub_atac_proj, 
                    trajectory = trajectory_name, 
                    colorBy= "cellColData", 
                    name = trajectory_name, 
                    embedding = "customUMAP", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-UMAP.pdf"), w = 4, h = 4)
p[[1]]
dev.off()

# Create heatmaps for each relevant matrix along the trajectory
# Get trajectories for each matrix
trajMM  <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL, trajectoryLabel = "FineNamedClust")
trajPM <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "PeakMatrix", trajectoryLabel = "FineNamedClust")
trajGS <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "GeneScoreMatrix", trajectoryLabel = "FineNamedClust")
trajGE <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "GeneExpressionMatrix", trajectoryLabel = "FineNamedClust")

# # Plot trajectories
genes.to.highlight <- c("FBXL7","NFIB","APLNR", "MEIS2", "VWF", "EDN1", "NFATC2", "SMAD1", "SERPINB1", 
                        "AQP1", "EDNRB", "LYVE1", "CHL1", 
                        "ITGA1", "S100A3", "CDH8", "ITGA4")

trajGE_genes <- gsub(".*:\\s*","",rownames(trajGE))

genes.idx <- which(trajGE_genes %in% genes.to.highlight)
trajGE_genes.to.highlight <- rownames(trajGE)[genes.idx]

p1 <- plotTrajectoryHeatmap(trajPM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), 
                            labelTop = 100)
p2 <- plotTrajectoryHeatmap(trajGS, 
                            pal = paletteContinuous(set = "horizonExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), 
                            labelTop = 100)
p3 <- plotTrajectoryHeatmap(trajGE, 
                            pal = paletteContinuous(set = "blueYellow"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), 
                            labelTop = 0,
                            labelMarkers = trajGE_genes.to.highlight
                            )
p4 <- plotTrajectoryHeatmap(trajMM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), 
                            labelTop = 100)

pdf(paste0(trajDir, "/TrajHeatmaps_", trajectory_name, ".pdf"), w = 8, h = 8)
p1
p2
p3
p4
dev.off()

# Endothelial trajectory (egCap -> eAero -> lAero)
# Define the endothelial trajectory based on optimal transport results
trajectory <- c("egCap", "eAero", "lAero")
trajectory_name <- "egCap_eAero_lAero_slingshot"

sub_atac_proj <- addSlingShotTrajectories(ArchRProj = sub_atac_proj, 
                               name = trajectory_name,
                               useGroups = trajectory, 
                               principalGroup = trajectory[1],
                               groupBy = "FineNamedClust", 
                               reducedDims = "PCA",
                               embedding = "customUMAP"
)

p <- plotTrajectory(sub_atac_proj, 
                    trajectory = paste0(trajectory_name, ".Curve1"), 
                    colorBy= "cellColData", 
                    name = paste0(trajectory_name, ".Curve1"), 
                    embedding = "customUMAP", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-UMAP.pdf"), w = 4, h = 4)
p[[1]]
dev.off()

# Create heatmaps for each relevant matrix along the trajectory
# Get trajectories for each matrix
trajMM  <- getTrajectory(ArchRProj = sub_atac_proj, name = paste0(trajectory_name, ".Curve1"), useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL, trajectoryLabel = "FineNamedClust")
trajPM <- getTrajectory(ArchRProj = sub_atac_proj, name = paste0(trajectory_name, ".Curve1"), useMatrix = "PeakMatrix", trajectoryLabel = "FineNamedClust")
trajGS <- getTrajectory(ArchRProj = sub_atac_proj, name = paste0(trajectory_name, ".Curve1"), useMatrix = "GeneScoreMatrix", trajectoryLabel = "FineNamedClust")
trajGE <- getTrajectory(ArchRProj = sub_atac_proj, name = paste0(trajectory_name, ".Curve1"), useMatrix = "GeneExpressionMatrix", trajectoryLabel = "FineNamedClust")

# # Plot trajectories
genes.to.highlight <- c("FBXL7","NFIB","APLNR", "MEIS2", "VWF", "EDN1", "NFATC2", "SMAD1", "SERPINB1", 
                        "AQP1", "EDNRB", "LYVE1",
                        "ITGA1", "S100A3", "CDH8", "ITGA4")

trajGE_genes <- gsub(".*:\\s*","",rownames(trajGE))

genes.idx <- which(trajGE_genes %in% genes.to.highlight)
trajGE_genes.to.highlight <- rownames(trajGE)[genes.idx]

p1 <- plotTrajectoryHeatmap(trajPM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), 
                            labelTop = 20)
p2 <- plotTrajectoryHeatmap(trajGS, 
                            pal = paletteContinuous(set = "horizonExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), 
                            labelTop = 20)
p3 <- plotTrajectoryHeatmap(trajGE, 
                            pal = paletteContinuous(set = "blueYellow"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), 
                            labelTop = 20,
                            labelMarkers = trajGE_genes.to.highlight
)
p4 <- plotTrajectoryHeatmap(trajMM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), 
                            labelTop = 20)

pdf(paste0(trajDir, "/TrajHeatmaps_", trajectory_name, ".pdf"), w = 8, h = 8)
p1
p2
p3
p4
dev.off()

# Correlated TF expression and motif activity along trajectory
# Correlate trajectories
corTraj_GE_MM <- correlateTrajectories(trajGE, trajMM)

trajGE2 <- trajGE[corTraj_GE_MM[["correlatedMappings"]]$name1,]
trajMM2 <- trajMM[corTraj_GE_MM[["correlatedMappings"]]$name2,]

trajCombined <- trajGE2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGE2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGE2))

ht1 <- plotTrajectoryHeatmap(trajGE2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, force = TRUE)

pdf(paste0(trajDir, "/corTraj-GE-MM-", trajectory_name,".pdf"))
ComplexHeatmap::draw(ht1 + ht2)
dev.off()

# Endothelial trajectory (egCap -> eArtr -> lArtr)
# Define the endothelial trajectory based on optimal transport results
trajectory <- c("egCap", "eArtr", "lArtr")
trajectory_name <- "egCap_eArtr_lArtr"

sub_atac_proj <- addTrajectory(ArchRProj = sub_atac_proj, 
                               name = trajectory_name, 
                               groupBy = "FineNamedClust", 
                               trajectory = trajectory, 
                               reducedDims = "PCA",
                               force = TRUE
                               #embedding = "customUMAP"
)

p <- plotTrajectory(sub_atac_proj, 
                    trajectory = trajectory_name, 
                    colorBy= "cellColData", 
                    name = trajectory_name, 
                    embedding = "customUMAP", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-UMAP.pdf"), w = 4, h = 4)
p[[1]]
dev.off()

# Create heatmaps for each relevant matrix along the trajectory
# Get trajectories for each matrix
trajMM  <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL, trajectoryLabel = "FineNamedClust")
trajPM <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "PeakMatrix", trajectoryLabel = "FineNamedClust")
trajGS <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "GeneScoreMatrix", trajectoryLabel = "FineNamedClust")
trajGE <- getTrajectory(ArchRProj = sub_atac_proj, name = trajectory_name, useMatrix = "GeneExpressionMatrix", trajectoryLabel = "FineNamedClust")

# # Plot trajectories
genes.to.highlight <- c("SOX6", "SOX5", "GJA5", "VEGFC", "DKK2", "FBLN5", "CA4", "MEF2C")

trajGE_genes <- gsub(".*:\\s*","",rownames(trajGE))

genes.idx <- which(trajGE_genes %in% genes.to.highlight)
trajGE_genes.to.highlight <- rownames(trajGE)[genes.idx]

p1 <- plotTrajectoryHeatmap(trajPM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), 
                            labelTop = 100)
p2 <- plotTrajectoryHeatmap(trajGS, 
                            pal = paletteContinuous(set = "horizonExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), 
                            labelTop = 100)
p3 <- plotTrajectoryHeatmap(trajGE, 
                            pal = paletteContinuous(set = "blueYellow"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), 
                            labelTop = 20#, 
                            #labelMarkers = trajGE_genes.to.highlight
)
p4 <- plotTrajectoryHeatmap(trajMM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), 
                            labelTop = 100)

pdf(paste0(trajDir, "/TrajHeatmaps_", trajectory_name, ".pdf"), w = 8, h = 8)
p1
p2
p3
p4
dev.off()

# Save ArchR SummarizedExperiment object for each trajectories
saveRDS(trajMM, paste0(trajDir, "/TrajMM_", trajectory_name, ".rds"))
saveRDS(trajPM, paste0(trajDir, "/TrajPM_", trajectory_name, ".rds"))
saveRDS(trajGS, paste0(trajDir, "/TrajGS_", trajectory_name, ".rds"))
saveRDS(trajGE, paste0(trajDir, "/TrajGE_", trajectory_name, ".rds"))

#trajMM <- readRDS(paste0(trajDir, "/TrajMM_", trajectory_name, ".rds"))
#trajPM <- readRDS(paste0(trajDir, "/TrajPM_", trajectory_name, ".rds"))
#trajGS <- readRDS(paste0(trajDir, "/TrajGS_", trajectory_name, ".rds"))
#trajGE <- readRDS(paste0(trajDir, "/TrajGE_", trajectory_name, ".rds"))

# save trajectory matrices used for plotting
trajPM_matrix <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), returnMatrix = T)
trajGS_matrix <- plotTrajectoryHeatmap(trajGS, pal = paletteContinuous(set = "horizonExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), returnMatrix = T)
trajGE_matrix <- plotTrajectoryHeatmap(trajGE, pal = paletteContinuous(set = "blueYellow"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), returnMatrix = T)
trajMM_matrix <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), returnMatrix = T)

saveRDS(trajPM_matrix, paste0(trajDir, "/TrajMM_", trajectory_name, "_matrix.rds"))
saveRDS(trajGS_matrix, paste0(trajDir, "/TrajPM_", trajectory_name, "_matrix.rds"))
saveRDS(trajGE_matrix, paste0(trajDir, "/TrajGS_", trajectory_name, "_matrix.rds"))
saveRDS(trajMM_matrix, paste0(trajDir, "/TrajGE_", trajectory_name, "_matrix.rds"))

# Positive TF regulators
# Identify deviant TF motifs
seGroupMotif <- getGroupSE(ArchRProj = sub_atac_proj, useMatrix = "MotifMatrix", groupBy = "FineNamedClust")
seGroupMotif <- seGroupMotif[,which(colnames(seGroupMotif) %in% trajectory)]
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

# Correlate motifmatrix and gene expression matrix
corGE_MM <- correlateMatrices(
  ArchRProj = sub_atac_proj,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PCA"
)

# Annotate each motif with the maximum delta observed between clusters from the initial deviant TF motifs
corGE_MM$maxDelta <- rowData(seZ)[match(corGE_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

# Identify positive TF regulators
# "we consider positive regulators as those TFs whose correlation between motif and gene score 
# (or gene expression) is greater than 0.5 with an adjusted p-value less than 0.01 and a maximum 
# inter-cluster difference in deviation z-score that is in the top quartile (Max TF Motif Delta)."
corGE_MM <- corGE_MM[order(abs(corGE_MM$cor), decreasing = TRUE), ]
corGE_MM <- corGE_MM[which(!duplicated(gsub("\\-.*","",corGE_MM[,"MotifMatrix_name"]))), ]
corGE_MM$TFRegulator <- "NO"
corGE_MM$TFRegulator[which(abs(corGE_MM$cor) > 0.5 & corGE_MM$padj < 0.01 & corGE_MM$maxDelta > quantile(corGE_MM$maxDelta, 0.75))] <- "YES"
sort(corGE_MM[corGE_MM$TFRegulator=="YES",1])

# Plot the positive TF regulators
p <- ggplot(data.frame(corGE_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  ggrepel::geom_label_repel(
    data = data.frame(corGE_MM[corGE_MM$TFRegulator=="YES",]), aes(x = cor, y = maxDelta, label = GeneExpressionMatrix_matchName), 
    size = 1.5,
    #nudge_x = 2,
    color = "black") +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGE_MM$maxDelta)*1.05)
  )

pdf(paste0(trajDir, "/PositiveTF_Regulators_", trajectory_name, ".pdf"), width = 5, height = 5)
p
dev.off()

# Correlated TF expression and motif activity along trajectory
# Correlate trajectories
corTraj_GE_MM <- correlateTrajectories(trajGE, trajMM)

trajGE2 <- trajGE[corTraj_GE_MM[["correlatedMappings"]]$name1,]
trajMM2 <- trajMM[corTraj_GE_MM[["correlatedMappings"]]$name2,]

trajCombined <- trajGE2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGE2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGE2))

ht1 <- plotTrajectoryHeatmap(trajGE2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, force = TRUE)

pdf(paste0(trajDir, "/corTraj-GE-MM-", trajectory_name,".pdf"))
ComplexHeatmap::draw(ht1 + ht2)
dev.off()

##########################################################################################
# Browser tracks and violin plots for Peak to gene linkages of trajectory correlated genes
##########################################################################################
driver.dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/3_trajectories/Endothelial/driver_genes"

clustOrder <- fineOrder %>% unlist() %>% unname()
celltypes <- c("lgCap", "lAero", "lArtr", "Veno")

driver.list <- list()
for (ct in celltypes) {
  csv <- readr::read_csv(paste0(driver.dir, "/", ct, "_driver_genes.csv"),
                         col_names = c("gene", "corr", "pval", "qval", "ci_low", "ci_high"), skip = 1)
  csv.filtered <- csv %>% drop_na()
  driver.list[[ct]] <- csv.filtered
}

driver.list.filtered <- lapply(driver.list, function(x){
  return(top_n(x, n = 50, wt = corr))
})

for (ct in celltypes) {
  genes <- driver.list.filtered[[ct]]$gene
  p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "FineNamedClust", 
    useGroups = clustsToPlot,
    pal = FineNamedClustCmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(7, 0.2, 1.25, 2.5),
    geneSymbol = genes, 
    loops = getPeak2GeneLinks(atac_proj),
    tileSize=500
  )
  plotPDF(plotList = p, 
          name = paste0("Endothelial_LineageDrivers_P2G_Tracks_", ct,".pdf"), 
          ArchRProj = atac_proj, 
          addDOC = FALSE, 
          width = 6, height = 6)
}

# Plot P2G tracks for specific lineages
clustsToPlot <- clustOrder[which(clustOrder %in% unique(sub_atac_proj$FineNamedClust))]
p <- plotBrowserTrack(
  ArchRProj = sub_atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustsToPlot,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = c("ITGA4", "FOXF1", "GJA5", "SOX6", "SOX5", "SOX13", "SOX17", "KLF12"), 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500
)

plotPDF(plotList = p,
        name = "BrowserTrackPlot-P2G-CandidateGenes.pdf",
        ArchRProj = atac_proj,
        addDOC = FALSE,
        width = 6, height = 6)

pdf(paste0(plotDir, "/BrowserTrackPlot-P2G-ITGA4.pdf"), w = 6, h = 6)
print(grid::grid.draw(p$ITGA4))
dev.off()

pdf(paste0(plotDir, "/BrowserTrackPlot-P2G-FOXF1.pdf"), w = 6, h = 6)
print(grid::grid.draw(p$FOXF1))
dev.off()

pdf(paste0(plotDir, "/BrowserTrackPlot-P2G-GJA5.pdf"), w = 6, h = 6)
print(grid::grid.draw(p$GJA5))
dev.off()

p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustsToPlot,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = c("LHX6", "HOXA5", "HMGA2", "DNMT1", "TEAD1", "CREB5", "NR2F2", "FOS", "JUN", "ATF4", "BCL6B", "FOXF1", "STAT4", "FOXP2", "ZNF331", "KLF9", "JUNB", "REL", "NFKB1", "MEF2A", "AR", "RELB", "ESM1"), 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500
)

plotPDF(plotList = p,
        name = "Endothelial_LineageDrivers_P2G_Tracks_CorrelatedTFs.pdf",
        ArchRProj = atac_proj,
        addDOC = FALSE,
        width = 5, height = 5)


# Violin plots of RNA expression for select genes

# WARNING: this may not work if you have already assigned the new p2glinks to the project
#label_genes <- c("ITGA4", "GJA5", "FOXF1", "FOXP1", "FOXP2")
label_genes <- c("GJA5", "ELN", "ITGA4", "ITGA3", "CHL1", "SOX5", "SOX6", "SOX13", "SOX17", "REL", "FOXF1", "FOXP2", "MEF2A")
GEmat <- getMatrixFromProject(atac_proj, useMatrix="GeneExpressionMatrix")
data_mat <- assays(GEmat)[[1]]
rownames(data_mat) <- rowData(GEmat)$name
sub_mat <- data_mat[label_genes,]

# These DO NOT match the order of the above matrix by default
grouping_data <- data.frame(cluster=factor(atac_proj$FineNamedClust, 
                                           ordered=TRUE))
rownames(grouping_data) <- getCellNames(atac_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

clustsToPlot <- clustOrder[which(clustOrder %in% c("egCap", "eAero", "lAero"))]

pList <- list()
for(gn in label_genes){
  gn <- "CHL1"
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% clustsToPlot,]
  df$cluster <- factor(df$cluster, levels = clustsToPlot)
  
  covarLabel <- "cluster"  
  #df <- filter(df, gene > 0)
  #df$gene <- log2(df$gene + 1)
  
  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=log(gene+1), fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    #+ geom_boxplot(aes(fill=cluster), adjust = 1.0, scale='width', na.rm = T, position=dodge)
    #+ geom_jitter(aes(group=Sample), size=0.025, 
    #  position=position_jitterdodge(seed=1, jitter.width=0.05, jitter.height=0.0, dodge.width=dodge_width))
    #+ stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
    # width=0.75, position=dodge,show.legend = FALSE)
    + scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=FineNamedClustCmap)
    + guides(fill=guide_legend(title=covarLabel), 
             colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  p
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Violin_Expression_byClust_P2G_Browser_Aero.pdf"), width=3, height=2)
pList
dev.off()

clustsToPlot <- clustOrder[which(clustOrder %in% c("egCap", "eArtr", "lArtr"))]

pList <- list()
for(gn in label_genes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% clustsToPlot,]
  df$cluster <- factor(df$cluster, levels = clustsToPlot)
  
  covarLabel <- "cluster"  
  
  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    #+ geom_jitter(aes(group=Sample), size=0.025, 
    #  position=position_jitterdodge(seed=1, jitter.width=0.05, jitter.height=0.0, dodge.width=dodge_width))
    #+ stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
    # width=0.75, position=dodge,show.legend = FALSE)
    + scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=FineNamedClustCmap)
    + guides(fill=guide_legend(title=covarLabel), 
             colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Violin_Expression_byClust_P2G_Browser_Artr.pdf"), width=3, height=2)
pList
dev.off()

##########################################################################################
# Linked Peaks for contribution score calculations
##########################################################################################
library(rtracklayer)

# Block scientific notations
options(scipen = 999)

# Peak directory 
peakDir <- paste0(wd, "/endothelial_p2g_peaks")
dir.create(peakDir, showWarnings = FALSE, recursive = TRUE)

genes <- c("ITGA4", "GJA5", "FOXF1", "SOX13", "SOX17", "EDNRB", "ITGA1", 
           "S100A3", "NRG3", "ITGA3", "CHL1", "IRS1", "SORBS1", "CD83",
           "SULF1", "ELN", "MMP16")

# Get peak set
peaks_GR <- getPeakSet(atac_proj)
genes_GR <- atac_proj@geneAnnotation$genes

# Get genes
filtered_genes_GR <- genes_GR[mcols(genes_GR)$symbol %in% genes]

# extend 50kb on both sides of genes
extended_genes_gr <- GRanges(
  seqnames = seqnames(filtered_genes_GR),
  ranges = IRanges(
    start = pmax(start(filtered_genes_GR) - 50000, 1),  # Ensure start doesn't go below 1
    end = end(filtered_genes_GR) + 50000
  ),
  strand = strand(filtered_genes_GR),
  mcols = mcols(filtered_genes_GR)  # Preserve metadata
)

# Peaks that are within 50kb upstream and downstream of the genes
filtered_peaks_GR <- findOverlaps(peaks_GR, extended_genes_gr)
filtered_peaks_GR <- peaks_GR[queryHits(filtered_peaks_GR)]

# Get blacklist regions used by chrombpnet
blacklist <- import.bed("/oak/stanford/groups/wjg/skim/projects/LDA/resources/chrombpnet/blacklist.bed.gz")
# Extend blacklist by 1057bp on both sides
blacklist_ext <- resize(blacklist, width = width(blacklist) + (1057 * 2), fix = "center")

filtered_peaks_GR <- filtered_peaks_GR[!filtered_peaks_GR %over% blacklist_ext]
filtered_peaks_GR$peakName <- paste0(seqnames(filtered_peaks_GR), "_", start(filtered_peaks_GR), "_", end(filtered_peaks_GR))

bed <- data.frame(
      seqnames = as.character(seqnames(filtered_peaks_GR)),
      starts = start(filtered_peaks_GR) - 251, #BED files are 0-indexed but granges are not
      ends = end(filtered_peaks_GR) + 250, # Convert to integer to prevent bed files being written in sci notation
      names = filtered_peaks_GR$peakName,
      scores = 0,
      strands = ".",
      col_6 = ".",
      col_7 = ".",
      col_8 = ".",
      summit = 500 #chrombpnet uses the 9th column as the summit position to calculate gc content. 500 since ArchR creates fixed width peaks that are 501bp long
      # and chrombpnet takese in 1000bp peak regions based on Laksshman paper
      )
write.table(bed, file = paste0(peakDir, "/endo_lineage_gene_peaks.bed"), quote = F, sep = "\t", row.names = F, col.names = F)

#For getting only P2G peaks
p2g_GR <- getPeak2GeneLinks(atac_proj, returnLoops = F)

peaks.list <- list()
for (gene in genes) {
  idxATAC <- dplyr::filter(as.data.frame(p2g_GR),
                           idxRNA == which(metadata(p2g_GR)$geneSet$name == gene)) %>%
    dplyr::select(idxATAC) %>%
    as.list() %>% unlist() %>% unname()
  peaks <- p2g_GR@metadata$peakSet[idxATAC]
  peaks.list[[gene]] <- peaks
}



##########################################################################################
# Sequence of promoter accessibility and gene expression
##########################################################################################

# Ensure trajectory has been added to sub_atac_proj
# Ensure driver.list has been loaded

# Aerocyte fate driver gene trendline
fate <- "lAero"
trajectory_name <- "egCap_eAero_lAero"
gene_set <- filter(driver.list[[fate]], corr > 0.5) %>% dplyr::select(gene) %>% as.list() %>% unname()

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "GeneExpressionMatrix",
                                name = paste0(fate, "_GE"),
                                features = gene_set)
p1 <- plotTrajectory(sub_atac_proj, 
                    trajectory = trajectory_name, 
                    colorBy= "cellColData", 
                    name = paste0(fate, "_GE1"), 
                    embedding = "customUMAP", addArrow = F)

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "GeneScoreMatrix",
                                name = paste0(fate, "_GS"),
                                features = gene_set)
p2 <- plotTrajectory(sub_atac_proj, 
                    trajectory = trajectory_name, 
                    colorBy= "cellColData", 
                    name = paste0(fate, "_GS1"), 
                    embedding = "customUMAP", addArrow = F)

peaks_gr <- getPeakSet(sub_atac_proj)
peaks_gr_near_gene <- peaks_gr[which(peaks_gr$nearestGene %in% unlist(gene_set)),]
peaks_gr_near_gene <- as.list(GRangesList(peaks_gr_near_gene))

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "PeakMatrix",
                                name = paste0(fate, "_PM"),
                                features = peaks_gr_near_gene)

p3 <- plotTrajectory(sub_atac_proj, 
                     trajectory = trajectory_name, 
                     colorBy= "cellColData", 
                     name = paste0(fate, "_PM1"), 
                     embedding = "customUMAP", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-", fate, "-fate_correlated_genes-trendline.pdf"), w = 4, h = 4)
p1[[2]]
p2[[2]]
p3[[2]]
dev.off()

# Arterial fate driver gene trendline
fate <- "lArtr"
trajectory_name <- "egCap_eArtr_lArtr"
gene_set <- filter(driver.list[[fate]], corr > 0.5) %>% dplyr::select(gene) %>% as.list() %>% unname()

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "GeneExpressionMatrix",
                                name = paste0(fate, "_GE"),
                                features = gene_set)
p1 <- plotTrajectory(sub_atac_proj, 
                     trajectory = trajectory_name, 
                     colorBy= "cellColData", 
                     name = paste0(fate, "_GE1"), 
                     embedding = "customUMAP", addArrow = F)

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "GeneScoreMatrix",
                                name = paste0(fate, "_GS"),
                                features = gene_set)
p2 <- plotTrajectory(sub_atac_proj, 
                     trajectory = trajectory_name, 
                     colorBy= "cellColData", 
                     name = paste0(fate, "_GS1"), 
                     embedding = "customUMAP", addArrow = F)

peaks_gr <- getPeakSet(sub_atac_proj)
peaks_gr_near_gene <- peaks_gr[which(peaks_gr$nearestGene %in% unlist(gene_set)),]
peaks_gr_near_gene <- as.list(GRangesList(peaks_gr_near_gene))

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "PeakMatrix",
                                name = paste0(fate, "_PM"),
                                features = peaks_gr_near_gene)

p3 <- plotTrajectory(sub_atac_proj, 
                     trajectory = trajectory_name, 
                     colorBy= "cellColData", 
                     name = paste0(fate, "_PM1"), 
                     embedding = "customUMAP", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-", fate, "-fate_correlated_genes-trendline.pdf"), w = 4, h = 4)
p1[[2]]
p2[[2]]
p3[[2]]
dev.off()

# GO term enrichment for 
source(paste0(scriptPath, "/GO_wrappers.R"))

fates <- c("lAero", "lArtr")
corr_threshold <- 0.5

count.mat <- Seurat::GetAssayData(object=rna_proj, slot="counts")
minUMIs <- 5
minCells <- 5
all_genes <- rownames(count.mat[rowSums(count.mat > minUMIs) > minCells,])

#all_genes <- c(driver.list$lAero$gene, driver.list$lArtr$gene) %>% unique() %>% sort()

gene_sets <- lapply(fates, function(f){
  driver.list[[f]] %>% dplyr::filter(corr > corr_threshold) %>% dplyr::select(gene) %>% unlist() %>% unname()
})

names(gene_sets) <- fates

GOresults <- lapply(fates, function(f){
  message(sprintf("Running GO enrichments on %s fate...", f))
  clust_genes <- gene_sets[[f]]
  upGO <- rbind(
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="BP") 
    #calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="CC")
  )
  upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
})

names(GOresults) <- paste0("fate_", fates)

# Plots of GO term enrichments:
pdf(paste0(trajDir, sprintf("/fate_driver_genes_GO_3termsBPonly_corr%s.pdf", corr_threshold)), width=10, height=2.5)
for(name in names(GOresults)){
  goRes <- GOresults[[name]]
  if(nrow(goRes)>1){
    print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
                       nterms=3, border_color="black", 
                       barwidth=0.85, title=name, barLimits=c(0, 15)))
  }
}
dev.off()


# Impute weights
sub_atac_proj <- addImputeWeights(sub_atac_proj, reducedDims = "PCA")

p <- plotTrajectory(sub_atac_proj,
                    colorBy = "GeneExpressionMatrix",
                    trajectory = trajectory_name,
                    name = "ITGA4",
                    continuousSet = "blueYellow",
                    embedding = "customUMAP")

p1 <- plotTrajectory(sub_atac_proj,
                     colorBy = "GeneScoreMatrix",
                     trajectory = trajectory_name,
                     name = "ITGA4",
                     continuousSet = "horizonExtra",
                     embedding = "customUMAP")

p2g_GR <- getPeak2GeneLinks(atac_proj, returnLoops = F)
idxATAC <- dplyr::filter(as.data.frame(p2g_GR), idxRNA == which(metadata(p2g_GR)$geneSet$name == "ITGA4")) %>% dplyr::select(idxATAC) %>% as.list() %>% unlist() %>% unname()
peaks <- lapply(idxATAC, function(idx){
  p2g_GR@metadata$peakSet[idx]
}) 
peaks <- unlist(as(peaks, "GRangesList"))
peaks <- as.list(GRangesList(peaks))

sub_atac_proj <- addModuleScore(sub_atac_proj,
                                useMatrix = "PeakMatrix",
                                name = "peaks_ITGA4",
                                features = peaks)

p2 <- plotTrajectory(sub_atac_proj,
                     colorBy = "cellColData",
                     trajectory = trajectory_name,
                     name = "peaks_ITGA41",
                     continuousSet = "solarExtra",
                     embedding = "customUMAP")

p2[[2]]

##########################################################################################
# Module scores of genes correlated with trajectory
##########################################################################################


DotPlot(sub_rna_proj, features = , group.by = "FineNamedClust")

intersect(sub(".*:", "", rownames(trajGE_matrix)), driver.list[["lAero"]]$gene)

sub_rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Endothelial/Endothelial.rds")

ct <- "lArtr"
set <- driver.list.filtered[[ct]]$gene
set <- list(c(set))

sub_rna_proj <- AddModuleScore(sub_rna_proj,
                               features = set,
                               name = "lArtr_traj")

FeaturePlot(sub_rna_proj, features = "lArtr_traj1")

#############################################################
# Identify regulatory targets of TFs 
#############################################################

# ChromVAR deviations matrix: (rows motif names x cols cell names)
motifMatrix <- getMatrixFromProject(sub_atac_proj, useMatrix="MotifMatrix")
deviationsMatrix <- assays(motifMatrix)$deviations

# GeneExpression Matrix: (rows gene names x cols cell names)
GEMatrix <- getMatrixFromProject(sub_atac_proj, useMatrix="GeneExpressionMatrix")
GEmat <- assays(GEMatrix)$GeneExpressionMatrix
rownames(GEmat) <- rowData(GEMatrix)$name
GEmat <- as(GEmat[Matrix::rowSums(GEmat) > 0,], "sparseMatrix") # Remove unexpressed genes

# Use only motifs that are 'TFRegulators' as determined by analysis above
GEMreg <- rownames(motifMatrix)[corGE_MM[corGE_MM$TFRegulator == "YES",]$MotifMatrix_idx]
regulators <- unique(GEMreg)

deviationsMatrix <- deviationsMatrix[regulators,]

# Identify pseudobulks for performing matrix correlations
knn_groups <- getLowOverlapAggregates(sub_atac_proj, target.agg=500, k=100, 
  overlapCutoff=0.8, dimReduc="PCA")

kgrps <- unique(knn_groups$group)

# GeneIntegrationMatrix
GEMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(GEmat[,use_cells])
  }) %>% do.call(cbind,.)
colnames(GEMatPsB) <- kgrps

# In rare instances, we can get pseudo-bulked genes that have zero averages
GEMatPsB <- GEMatPsB[Matrix::rowSums(GEMatPsB) > 0,]

# DeviationsMatrix
DevMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(deviationsMatrix[,use_cells])
  }) %>% do.call(cbind,.)
colnames(DevMatPsB) <- kgrps

# Perform chromVAR deviations toRNA correlation analysis:
start <- Sys.time()
geneCorMat <- cor2Matrices(DevMatPsB, GEMatPsB)
colnames(geneCorMat) <- c("motifName", "symbol", "Correlation", "FDR")
end <- Sys.time()
message(sprintf("Finished correlations in %s minutes.", round((end  - start)/60.0, 2)))

allGenes <- rownames(GEMatPsB) %>% sort() # Already filtered to only expressed genes

# Get locations of motifs of interest:
motifPositions <- getPositions(sub_atac_proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# Get peak to gene GR
p2gGR <- getP2G_GR(sub_atac_proj, corrCutoff = 0.45, filtNA = TRUE) #######################rerun this with correlation filter....


calculateLinkageScore <- function(motifLocs, p2gGR){
  # Calculate Linkage Score (LS) for each gene in p2gGR with regards to a motif location GR
  ###################################
  # For a given gene, the LS = sum(corr peak R2 * motifScore)
  ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
  olGenes <- p2gGR[to(ol)]
  olGenes$motifScore <- motifLocs[from(ol)]$score
  olGenes$R2 <- olGenes$Correlation**2 # All p2g links here are already filtered to only be positively correlated
  LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
  LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
  LSdf$rank <- 1:nrow(LSdf)
  return(LSdf)
}

calculateMotifEnrichment <- function(motifLocs, p2gGR){
  # Calculate Motif enrichment per gene
  ###################################
  # For a given gene, calculate the hypergeometric enrichment of motifs in 
  # linked peaks (generally will be underpowered)
  motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
  m <- length(motifP2G) # Number of possible successes in background
  n <- length(p2gGR) - m # Number of non-successes in background

  motifLinks <- motifP2G$symbol %>% getFreqs()
  allLinks <- p2gGR$symbol %>% getFreqs()
  df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
  df$motifLinks[is.na(df$motifLinks)] <- 0
  df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
  df <- df[order(df$mLog10pval, decreasing=TRUE),]
  df$symbol <- rownames(df)
  return(df)
}

# plot all TF regulators
regPlotDir <- paste0(plotDir, "/TFregulatorPlots")
dir.create(regPlotDir, showWarnings = FALSE, recursive = TRUE)

markerGenes  <- c(
  featureSets
) %>% unlist() %>% unname() %>% unique()

# Get list of genes we want to highlight (e.g. genes involved in HF development)
# library(org.Hs.eg.db)
# library(GO.db)
# go_id = GOID(GOTERM[Term(GOTERM) == "cornification"])
# allegs = get(go_id, org.Hs.egGO2ALLEGS)
# cornification_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()
# 
# go_id = GOID(GOTERM[Term(GOTERM) == "hemidesmosome assembly"])
# allegs = get(go_id, org.Hs.egGO2ALLEGS)
# hemidesmosome_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()
# 
# go_id = GOID(GOTERM[Term(GOTERM) == "hair follicle development"])
# allegs = get(go_id, org.Hs.egGO2ALLEGS)
# hfdev_genes = mget(allegs,org.Hs.egSYMBOL) %>% unlist() %>% unname() %>% unique() %>% sort()
# 
# markerGenes <- c(cornification_genes, hemidesmosome_genes, hfdev_genes, markerGenes) %>% unique() %>% sort()

# Store results for each TF
res_list <- list()

###########################################
for(motif in regulators){
  motif_short <- strsplit(motif,"_")[[1]][1]
  # First get motif positions
  motifLocs <- motifGR[motifGR$motifName == motif]
  # Calculate Linkage Score for motif
  LS <- calculateLinkageScore(motifLocs, p2gGR)
  # Get just genes correlated to motif
  motifGeneCorDF <- geneCorMat[geneCorMat$motifName == motif,]
  plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
  # Calculate motif enrichment per gene
  ME <- calculateMotifEnrichment(motifLocs, p2gGR)
  plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
  plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
  plot_df$toLabel <- "NO"
  topN <- 5
  plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
  plot_df$rank_LS <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
  plot_df$rank_Corr <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
  plot_df$rank_Pval <- 1:nrow(plot_df)
  plot_df$toLabel[1:10] <- "YES"
  plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
  plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
  plot_df$toLabel[1:topN] <- "YES"
  # Label any marker genes in window of interest
  LS_window <- quantile(plot_df$LS, 0.8)
  corr_window <- 0.25
  pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
  neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
  if(nrow(plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]) > 0){
    plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
  }
  res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
  # Save dataframe of results
  save_df <- plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes),c(1:5)]
  save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
  saveRDS(save_df, paste0(regPlotDir, sprintf("/regulatory_targets_%s.rds", motif_short)))
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
  # Label motif as well
  plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
  plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
  # Threshold pvalue for plotting
  maxPval <- 5
  plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
  #Plot results
  p <- (
    ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
      #+ geom_point(size = 2)
      + ggrastr::geom_point_rast(size=2)
      + ggrepel::geom_text_repel(
          data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
          #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
          size=2,
          point.padding=0, # additional pading around each point
          box.padding=0.5,
          min.segment.length=0, # draw all line segments
          max.overlaps=Inf, # draw all labels
          #nudge_x = 2,
          color="black"
      ) 
      + geom_vline(xintercept=0, lty="dashed") 
      + geom_vline(xintercept=corr_window, lty="dashed", color="red")
      + geom_vline(xintercept=-corr_window, lty="dashed", color="red")
      + geom_hline(yintercept=LS_window, lty="dashed", color="red")
      + theme_BOR(border=FALSE)
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio=1.0,
              #legend.position = "none", # Remove legend
              axis.text.x = element_text(angle=90, hjust=1))
      + ylab("Linkage Score") 
      + xlab("Motif Correlation to Gene") 
      + scale_color_gradientn(colors=cmaps_BOR$zissou, limits=c(0, maxPval))
      + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
      + ggtitle(sprintf("%s putative targets", motif_short))
      )
  # Positively regulated genes:
  upGO <- rbind(
    calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="MF")
    )
  upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  up_go_plot <- topGObarPlot(upGO, cmap=cmaps_BOR$comet, nterms=6, border_color="black", 
    barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(pos_top_genes)), enrichLimits=c(0, 6))
  # Negatively regulated genes:
  downGO <- rbind(
    calcTopGo(allGenes, interestingGenes=neg_top_genes, nodeSize=5, ontology="BP"), 
    calcTopGo(allGenes, interestingGenes=neg_top_genes, nodeSize=5, ontology="MF")
    )
  downGO <- downGO[order(as.numeric(downGO$pvalue), decreasing=FALSE),]
  down_go_plot <- topGObarPlot(downGO, cmap=cmaps_BOR$comet, nterms=6, border_color="black", 
    barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(neg_top_genes)), enrichLimits=c(0, 6))
  pdf(paste0(regPlotDir, sprintf("/%s_LS.pdf", motif_short)), width=8, height=6)
  print(p)
  print(up_go_plot)
  print(down_go_plot)
  dev.off()
}
###########################################


##########################################################################################
# Peak2Gene Link Heatmap
##########################################################################################
p2gGR <- getP2G_GR(sub_atac_proj)
p2gFreqs <- getFreqs(p2gGR$symbol)

x <- 1:length(p2gFreqs)
rank_df <- data.frame(npeaks=p2gFreqs, rank=x)

p <- plotBrowserTrack(
  ArchRProj = sub_atac_proj,
  groupBy = "FineNamedClust",
  geneSymbol = c("EDNRB", "GJA5", "DKK2", "ACKR1", "SOX5"),
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(sub_atac_proj)
)

grid::grid.draw(p$SOX5)

p <- plotPeak2GeneHeatmap(
  ArchRProj = sub_atac_proj, 
  corCutOff = 0.5,
  groupBy = "FineNamedClust", 
  k = 4)

pdf(paste0(plotDir, "/P2G_Heatmap_Endothelial.pdf"))
p
dev.off()

p2g_df <- plotPeak2GeneHeatmap(
  ArchRProj = sub_atac_proj, 
  corCutOff = 0.5,
  groupBy = "FineNamedClust", 
  k = 4, returnMatrices = T)

saveRDS(p2g_df, paste0(plotDir, "/P2G_Heatmap_Endothelial_Kmeans_k4_df.rds"))

##########################################################################################
# Motifs around target TFs
##########################################################################################
# TRial
tfs <- c("NR2F1", "TWIST2", "GATA6", "ETS1", "STAT4", "CREB5", "NFATC2", "PRRX2", "LEF1" ,"MEF2C", "MEF2B", "EBF2")

# Get locations of motifs of interest:
motifPositions <- getPositions(sub_atac_proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# First get motif positions
#motifLocs <- motifGR[motifGR$motifName %in% tfs]

# Get peak2gene linkages for the target TF
p2gGR <- getP2G_GR(sub_atac_proj, corrCutoff = 0.45, filtNA = TRUE)
p2gGR_metadata <- mcols(p2gGR)

motif <- "FOXP2"
p2g_motif <- p2gGR[which(p2gGR_metadata$symbol %in% motif),]

# Which motifs are present in the peak2gene linkages of the target TF
ol <- findOverlaps(motifGR, p2g_motif, maxgap = 0 , type = c("any"), ignore.strand = TRUE)

motifGR[from(ol)]$motifName %>% as.vector()

##########################################################################################
# Violin plots for Expression of genes for PROX1 browser tracks
##########################################################################################

label_genes <- c("NFATC1", "HNF1A" ,"HNF1B", "GATA2", "FOXC2")
GEmat <- getMatrixFromProject(atac_proj, useMatrix="GeneExpressionMatrix")
data_mat <- assays(GEmat)[[1]]
rownames(data_mat) <- rowData(GEmat)$name
sub_mat <- data_mat[label_genes,]

# These DO NOT match the order of the above matrix by default
grouping_data <- data.frame(cluster=factor(atac_proj$FineNamedClust, 
                                           ordered=TRUE))
rownames(grouping_data) <- getCellNames(atac_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

clustsToPlot <- clustOrder[which(clustOrder %in% c("Lymp", "PNEC1", "PNEC2", "PNEC3"))]

pList <- list()
for(gn in label_genes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% clustsToPlot,]
  df$cluster <- factor(df$cluster, levels = clustsToPlot)
  
  covarLabel <- "cluster"  
  #df <- filter(df, gene > 0)
  #df$gene <- log2(df$gene + 1)
  
  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=log(gene+1), fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    #+ geom_boxplot(aes(fill=cluster), adjust = 1.0, scale='width', na.rm = T, position=dodge)
    #+ geom_jitter(aes(group=Sample), size=0.025, 
    #  position=position_jitterdodge(seed=1, jitter.width=0.05, jitter.height=0.0, dodge.width=dodge_width))
    #+ stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
    # width=0.75, position=dodge,show.legend = FALSE)
    + scale_color_manual(values=FineNamedClustCmap, limits=names(FineNamedClustCmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=FineNamedClustCmap)
    + guides(fill=guide_legend(title=covarLabel), 
             colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  p
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Violin_Expression_Prox1BrowserTrack.pdf"), width=3, height=2)
pList
dev.off()
