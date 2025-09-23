#!/usr/bin/env Rscript

# Figure 4 panels ####

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
source(paste0(scriptPath, "/sample_metadata.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/Figure04_Epithelial"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

# set plot directory
plotDir <- paste0(wd)

# Misc options
addArchRGenome("hg38")
pointSize <- 1

# Preparing Data #####
# Load data objects
atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_v2", force = TRUE)
rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# Color Maps
compartmentCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/lungClusterColors.rds") %>% unlist()
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds") %>% unlist()
sample_cmap <- readRDS(paste0("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds"))
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")

# Load subprojects
sub_rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Epithelial/Epithelial.rds")
sub_atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Epithelial")

Idents(sub_rna_proj) <- "FineNamedClust"


# UMAPs ##########

# Add full RNA object UMAP to subclustered object
umap <- rna_proj@reductions$umap
sub_rna_proj@reductions$umap_full <- umap

# Add PCs from full RNA object to subclustered object
pca <- rna_proj@reductions$pca
idx <- which(rownames(pca) %in% colnames(sub_rna_proj))
pca@cell.embeddings <- pca@cell.embeddings[idx,]
sub_rna_proj@reductions$pca <- pca

# UMAP on the subsetted object
# UMAP parameters
umapDims <- 40
umapNeighbors <- 50
umapMinDist <- 0.5
umapDistMetric <- "cosine"

sub_rna_proj <- RunUMAP(object = sub_rna_proj, 
        reduction = "pca",
        assay = "RNA",
        dims = 1:umapDims,
        n.neighbors = umapNeighbors,
        min.dist = umapMinDist, force = TRUE)

# Select straggler cells to remove that are not on the UMAP
plot <- DimPlot(sub_rna_proj, reduction = "umap_full", group.by = "FineNamedClust")
HoverLocator(plot = plot, information = FetchData(sub_rna_proj, vars = c("ident", "FineNamedClust")))
selected.cells <- CellSelector(plot = plot)
selected.cells.1 <- CellSelector(plot = plot)
selected.cells.2 <- CellSelector(plot = plot)

cellsToRemove <- c(selected.cells, selected.cells.1, selected.cells.2, colnames(sub_rna_proj[,which(sub_rna_proj$FineNamedClust == "Schw")]))
saveRDS(cellsToRemove, file=paste0(plotDir, "/cellsToRemove_epithelial_subcluster.rds"))
cellsToRemove <- readRDS(file=paste0(plotDir, "/cellsToRemove_epithelial_subcluster.rds"))

filt_sub_rna_proj <- sub_rna_proj[,which(colnames(sub_rna_proj) %ni% cellsToRemove)]
filt_sub_rna_proj@reductions$umap <- filt_sub_rna_proj@reductions$umap_full
saveRDS(filt_sub_rna_proj, file = "/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Epithelial/filtered_Epithelial.rds")
filt_sub_rna_proj <- readRDS(file = "/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Epithelial/filtered_Epithelial.rds")

keep_cbs <- colnames(filt_sub_rna_proj)
stringr::str_sub(keep_cbs, -19, -19) <- "#"
keep_cbs <- keep_cbs[which(keep_cbs %in% getCellNames(sub_atac_proj))]

sub_atac_proj <- addGroupCoverages(sub_atac_proj, groupBy = "FineNamedClust", force = TRUE)

# Temp for now do not write or save anything due to GroupCoverages directory not pointing to the correct directory
filt_sub_atac_proj <- sub_atac_proj[keep_cbs,]

# Add umap from filtered full RNA proj to sub ArchR Project
umap <- filt_sub_rna_proj@reductions$umap_full@cell.embeddings
CBs <- rownames(umap)
str_sub(CBs, -19, -19) <- "#"
rownames(umap) <- CBs
umap <- data.frame(umap)

# Only keep cells in the atac project 
umap <- umap %>% filter(row.names(umap) %in% getCellNames(filt_sub_atac_proj))
colnames(umap) <- c("custom#UMAP1", "custom#UMAP2")
# Add seurat umap coordinates into ArchR
filt_sub_atac_proj@embeddings$customUMAP <- SimpleList(df = umap, params = list())

# Add PCs from full RNA object to subclustered atac project
pca <- Embeddings(rna_proj, reduction = "pca")
# convert cell names for ArchR
str_sub(rownames(pca), -19, -19) <- "#"
# clean up pca matrix for ArchR Project
pca <- pca[which(rownames(pca) %in% getCellNames(filt_sub_atac_proj)),]
str_sub(rownames(pca), -19, -19) <- "#"

# add RNA PCs to ArchR Project
filt_sub_atac_proj@reducedDims[["PCA"]] <- SimpleList(
  matDR = pca,
  date = Sys.time(),
  assay.used = rna_proj@reductions$pca@assay.used,
  scaleDims = NA,
  corToDepth = NA
)

umapPlots <- list()
# FineNamedClust on RNA clustering:
umapDF <- buildUMAPdfFromArchR(sub_atac_proj, cellColData="FineNamedClust", embeddingName="customUMAP_full")
readr::write_tsv(umapDF, file = paste0(plotDir, "/FineNamedClust_on_RNA.tsv"))
umapPlots[["FineNamedClust_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=FineNamedClustCmap, 
                                                 namedColors=TRUE, point_size=pointSize, covarLabel="FineNamedClust_on_RNA", useRaster=TRUE)

# Age on RNA clustering:
umapDF <- buildUMAPdfFromArchR(sub_atac_proj, cellColData="age", embeddingName="customUMAP_full")
readr::write_tsv(umapDF, file = paste0(dir.name, "/Ages_on_RNA.tsv"))
umapPlots[["Ages_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=gest_age_cmap, 
                                           namedColors=TRUE, point_size=pointSize, covarLabel="Age_on_RNA", useRaster=TRUE)

plot.name <- "UMAP_RNA_FineNamedClust"
dir.name <- paste0(wd, "/", plot.name)
dir.create(dir.name, showWarnings = FALSE, recursive = TRUE)

pdf(paste0(plotDir,"/",plot.name,".pdf"), width=7, height=5)
umapPlots
dev.off()


##########################################################################################
# Differential expression testing between AT2l and AT1l
##########################################################################################
library(DESeq2)

DE_dir <- paste0(wd, "/DE_testing")
dir.create(DE_dir, showWarnings = FALSE, recursive = TRUE)

# pseudobulk the counts based on sample and cell type
pseudo_rna <- AggregateExpression(sub_rna_proj, 
                                  assays = "RNA", 
                                  return.seurat = T, 
                                  group.by = c("Sample", "FineNamedClust"), slot = "counts")
pseudo_rna$Sample <- gsub("_(.*)", "", Cells(pseudo_rna))
pseudo_rna$FineNamedClust <- gsub(".*_", "", Cells(pseudo_rna)) 
Idents(pseudo_rna) <- "FineNamedClust"

# DE using DESeq2 
pseudo_rna_sub <- subset(pseudo_rna, ident = c("AT2l", "AT1l"))
data.use <- GetAssayData(pseudo_rna_sub, "counts", "RNA")

# Add metadata information on sample (batch)
group.info <- data.frame(celltype = as.factor(gsub(".*_", "", unname(colnames(data.use)))),
                         samplename = unname(colnames(data.use)), 
                         sample = gsub("_(.*)", "", unname(colnames(data.use)))
)
group.info$celltype <- relevel(group.info$celltype, ref = "AT2l") # set reference cell type to use for fold change

dds1 <- DESeq2::DESeqDataSetFromMatrix(
  countData = data.use,
  colData = group.info,
  design = ~ celltype + sample
)

# Run DESeq2
dds1 <- DESeq(dds1)
resultsNames(dds1)
resLFC <- lfcShrink(dds1, coef="celltype_AT1l_vs_AT2l", type="apeglm")
resOrdered <- resLFC[order(resLFC$pvalue),]

saveRDS(dds1, paste0(DE_dir, "/DDS.rds"))
saveRDS(resOrdered, paste0(DE_dir, "/resOrdered.rds"))
resOrdered <- readRDS(paste0(DE_dir, "/resOrdered.rds"))

res.df <- resOrdered %>% as.data.frame()
res.df$gene <- rownames(res.df)
readr::write_csv(res.df, paste0(DE_dir, "/AT1lvsAT2l_DE_genes.csv"))
readr::write_tsv(res.df, paste0(DE_dir, "/AT1lvsAT2l_DE_genes.tsv"))

# Plot MA Plot of differential GE testing
l2fc_threshold <- 1.0
padj_threshold <- 0.01
geneToLabel <- c("AGER", "SFTPC", "WNT7A", "GPC5", "MYO16", "AQP5", "MYRF", "SFTPA1")
res.df$label <- NA
idxToLabel <- which(res.df$gene %in% geneToLabel)
res.df$label[idxToLabel] <- res.df$gene[idxToLabel]

res.df$color <- "grey78"
res.df$color[which(res.df$log2FoldChange > l2fc_threshold & res.df$padj < padj_threshold)] <- FineNamedClustCmap[["AT1l"]]
res.df$color[which(res.df$log2FoldChange < -l2fc_threshold & res.df$padj < padj_threshold)] <- FineNamedClustCmap[["AT2l"]]

res.plot <- res.df[,c("baseMean", "log2FoldChange", "label", "color")]

point_size = 1.5
color_col <- "color"

pdf(paste0(DE_dir, "/MA_plot_AT1lvsAT2l_l2fc", l2fc_threshold,"_padj", padj_threshold,".pdf"))
ggplot(data.frame(res.plot), aes(x=log2(res.plot[,1]), y=res.plot[,2])) +
  geom_point_rast(size=point_size, color = res.plot[,color_col]) +
  ggrepel::geom_label_repel(
    data = data.frame(res.plot[!is.na(res.plot$label),]), 
    aes(x = log2(baseMean), y = log2FoldChange, label = label), 
    size = 3, box.padding = 0.2,
    #nudge_x = 2,
    color = "black") +
  geom_hline(yintercept=c(-l2fc_threshold, l2fc_threshold), linetype="dashed") +
  xlab(paste0("log2",colnames(res.plot)[1])) +
  ylab(colnames(res.plot)[2]) +
  theme_BOR(border=TRUE) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        plot.margin=unit(c(0.25,1,0.25,1), "cm"), 
        aspect.ratio=1,
        #legend.position = "none", # Remove legend
        #axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  xlim(0, 15)
dev.off()

saveRDS(res.plot, file = paste0(DE_dir, "/MA_plot_AT1lvsAT2l_l2fc.rds"))

##########################################################################################
# Differential Peaks between AT2l and AT1l
##########################################################################################

diffPeaks <- getMarkerFeatures(
  ArchRProj = sub_atac_proj, 
  useMatrix = "PeakMatrix",
  groupBy = "FineNamedClust",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AT1l",
  bgdGroups = "AT2l"
)

pma <- plotMarkers(seMarker = diffPeaks, 
                   name = "AT1l", 
                   cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", 
                   plotAs = "MA")

pdf(paste0(plotDir, "/AT1lvsAT2l_diffPeaks_MAplot.pdf"))
pma
dev.off()

motifs_AT1l <- peakAnnoEnrichment(
  seMarker = diffPeaks,
  ArchRProj = sub_atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifs_AT2l <- peakAnnoEnrichment(
  seMarker = diffPeaks,
  ArchRProj = sub_atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= 0.5"
)

df1 <- data.frame(TF = rownames(motifs_AT1l), 
                 mlog10Padj = assay(motifs_AT1l)[,1],
                 celltype = "AT1l")
df1 <- df1[order(df1$mlog10Padj, decreasing = TRUE),]
df1$rank <- seq_len(nrow(df1))

df2 <- data.frame(TF = rownames(motifs_AT2l), 
                  mlog10Padj = assay(motifs_AT2l)[,1],
                  celltype = "AT2l")
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))

gg1 <- ggplot(df1, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df1[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

gg2 <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df2[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf(paste0(DE_dir, "/rankplot_motifenrich_diffPeaks_AT1l_AT2l.pdf"), w = 4, h = 4)
gg1 + ggtitle("AT1l")
gg2 + ggtitle("AT2l")
dev.off()


# Correlate motifmatrix and gene expression matrix ###
alveolar_atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Alveolar")

# Import Seurat PCA embeddings for RNA
pca <- Embeddings(rna_proj, reduction = "pca")
str_sub(rownames(pca), -19, -19) <- "#" # convert cell names for ArchR
pca <- pca[getCellNames(alveolar_atac_proj),]

# add RNA PCs to ArchR Project
alveolar_atac_proj@reducedDims[["PCA"]] <- SimpleList(
  matDR = pca,
  date = Sys.time(),
  assay.used = rna_proj@reductions$pca@assay.used,
  scaleDims = NA,
  corToDepth = NA
)

corGE_MM <- correlateMatrices(
  ArchRProj = alveolar_atac_proj,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PCA"
)

corGE_MM <- corGE_MM[order(abs(corGE_MM$cor), decreasing = TRUE), ]
corGE_MM <- corGE_MM[which(!duplicated(gsub("\\-.*","",corGE_MM[,"MotifMatrix_name"]))), ] # remove AS and other non-gene variants of TF

res.df <- readr::read_tsv(paste0(DE_dir, "/AT1lvsAT2l_DE_genes.tsv"))

l2fc_threshold <- 1
padj_threshold <- 0.01
plot.df <- left_join(as.data.frame(corGE_MM), as.data.frame(res.df), by = c("GeneExpressionMatrix_name" = "gene"))
plot.df$label <- "NO"
plot.df$label[which(plot.df$cor > 0.0 & plot.df$padj.x < padj_threshold & abs(plot.df$log2FoldChange) > l2fc_threshold & plot.df$padj.y < padj_threshold)] <- "YES"

pdf(paste0(DE_dir, "/AT_Motif_Activity_Expression_Correlation_By_AT1lvsAT2l_DiffExp.pdf"), w = 6, h = 6)
ggplot(data.frame(plot.df), aes(cor, log2FoldChange, color = label)) +
  geom_point() + 
  ggrepel::geom_label_repel(
    data = data.frame(plot.df[plot.df$label=="YES",]), 
    aes(x = cor, y = log2FoldChange, label = GeneExpressionMatrix_matchName), 
    size = 4,
    #nudge_x = 2,
    color = "black") +
  theme_ArchR() +
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation of TF Motif Accessibility to expression") +
  ylab("Log2FoldChange AT1l vs. AT2l") +
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(-3, 3), expand = c(0,0)) +
  geom_hline(yintercept = c(-l2fc_threshold, l2fc_threshold), linetype="dashed")
dev.off()

##########################################################################################
# Trajectory analyses
##########################################################################################

# Trajectory directory
trajDir <- paste0(wd, "/trajectories")
dir.create(trajDir, showWarnings = FALSE, recursive = TRUE)

# Ciliated cell trajectory #

# Define the trajectory based on optimal transport results
trajectory <- c("APr2", "eCili", "lCili")
trajectory_name <- paste(trajectory, collapse = "_")

filt_sub_atac_proj <- addTrajectory(ArchRProj = filt_sub_atac_proj, 
                               name = trajectory_name, 
                               groupBy = "FineNamedClust", 
                               trajectory = trajectory, 
                               reducedDims = "PCA",
                               force = TRUE,
                               embedding = "customUMAP_full"
)

p <- plotTrajectory(filt_sub_atac_proj, 
                    trajectory = trajectory_name, 
                    colorBy= "cellColData", 
                    name = trajectory_name, 
                    embedding = "customUMAP_full", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-UMAP.pdf"), w = 4, h = 4)
p[[1]]
dev.off()

# Create heatmaps for each relevant matrix along the trajectory
# Get trajectories for each matrix
trajMM  <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL, trajectoryLabel = "FineNamedClust")
trajPM <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "PeakMatrix", trajectoryLabel = "FineNamedClust")
trajGS <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "GeneScoreMatrix", trajectoryLabel = "FineNamedClust")
trajGE <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "GeneExpressionMatrix", trajectoryLabel = "FineNamedClust")

# # Plot trajectories
#genes.to.highlight <- c("SCGB3A2", "SOX2", "FOXJ1", "DNAH12", "DNAH9")
#trajGE_genes <- gsub(".*:\\s*","",rownames(trajGE))
#genes.idx <- which(trajGE_genes %in% genes.to.highlight)
#trajGE_genes.to.highlight <- rownames(trajGE)[genes.idx]

p1 <- plotTrajectoryHeatmap(trajPM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), 
                            labelTop = 50)
p2 <- plotTrajectoryHeatmap(trajGS, 
                            pal = paletteContinuous(set = "horizonExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), 
                            labelTop = 50)
p3 <- plotTrajectoryHeatmap(trajGE, 
                            pal = paletteContinuous(set = "blueYellow"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), 
                            labelTop = 25)
p4 <- plotTrajectoryHeatmap(trajMM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), 
                            labelTop = 50)

pdf(paste0(trajDir, "/TrajHeatmaps_", trajectory_name, ".pdf"), w = 8, h = 8)
p1
p2
p3
p4
dev.off()

# Correlated TF expression and motif activity along trajectory
# Correlate trajectories
corTraj_GE_MM <- correlateTrajectories(trajGE, trajMM, )
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

# Positive TF regulators for this trajectory
# Identify deviant TF motifs
seGroupMotif <- getGroupSE(ArchRProj = filt_sub_atac_proj, useMatrix = "MotifMatrix", groupBy = "FineNamedClust")
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

# Trajectory directory
trajDir <- paste0(wd, "/trajectories")
dir.create(trajDir, showWarnings = FALSE, recursive = TRUE)

# AT2l cell trajectory #######
# Define the trajectory based on optimal transport results
trajectory <- c("TiP1", "TiP2", "AT2l")
trajectory_name <- paste(trajectory, collapse = "_")

filt_sub_atac_proj <- addTrajectory(ArchRProj = filt_sub_atac_proj, 
                                    name = trajectory_name, 
                                    groupBy = "FineNamedClust", 
                                    trajectory = trajectory, 
                                    reducedDims = "PCA",
                                    force = TRUE,
                                    embedding = "customUMAP_full"
)

p <- plotTrajectory(filt_sub_atac_proj, 
                    trajectory = trajectory_name, 
                    colorBy= "cellColData", 
                    name = trajectory_name, 
                    embedding = "customUMAP_full", addArrow = F)

pdf(paste0(trajDir, "/Traj-", trajectory_name, "-UMAP.pdf"), w = 4, h = 4)
p[[1]]
dev.off()

# Create heatmaps for each relevant matrix along the trajectory
# Get trajectories for each matrix
trajMM  <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL, trajectoryLabel = "FineNamedClust")
trajPM <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "PeakMatrix", trajectoryLabel = "FineNamedClust")
trajGS <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "GeneScoreMatrix", trajectoryLabel = "FineNamedClust")
trajGE <- getTrajectory(ArchRProj = filt_sub_atac_proj, name = trajectory_name, useMatrix = "GeneExpressionMatrix", trajectoryLabel = "FineNamedClust")

# # Plot trajectories
#genes.to.highlight <- c("SCGB3A2", "SOX2", "FOXJ1", "DNAH12", "DNAH9")
#trajGE_genes <- gsub(".*:\\s*","",rownames(trajGE))
#genes.idx <- which(trajGE_genes %in% genes.to.highlight)
#trajGE_genes.to.highlight <- rownames(trajGE)[genes.idx]

p1 <- plotTrajectoryHeatmap(trajPM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), 
                            labelTop = 50)
p2 <- plotTrajectoryHeatmap(trajGS, 
                            pal = paletteContinuous(set = "horizonExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), 
                            labelTop = 50)
p3 <- plotTrajectoryHeatmap(trajGE, 
                            pal = paletteContinuous(set = "blueYellow"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), 
                            labelTop = 25)
p4 <- plotTrajectoryHeatmap(trajMM, 
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), 
                            labelTop = 50)

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

# save trajectory matrices used for plotting
trajPM_matrix <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)), returnMatrix = T)
trajGS_matrix <- plotTrajectoryHeatmap(trajGS, pal = paletteContinuous(set = "horizonExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGS)$label)), returnMatrix = T)
trajGE_matrix <- plotTrajectoryHeatmap(trajGE, pal = paletteContinuous(set = "blueYellow"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGE)$label)), returnMatrix = T)
trajMM_matrix <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)), returnMatrix = T)

saveRDS(trajPM_matrix, paste0(trajDir, "/TrajMM_", trajectory_name, "_matrix.rds"))
saveRDS(trajGS_matrix, paste0(trajDir, "/TrajPM_", trajectory_name, "_matrix.rds"))
saveRDS(trajGE_matrix, paste0(trajDir, "/TrajGS_", trajectory_name, "_matrix.rds"))
saveRDS(trajMM_matrix, paste0(trajDir, "/TrajGE_", trajectory_name, "_matrix.rds"))

# Correlated TF expression and motif activity along trajectory
# Correlate trajectories
corTraj_GE_MM <- correlateTrajectories(trajGE, trajMM, )
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

# Positive TF regulators for this trajectory
# Identify deviant TF motifs
seGroupMotif <- getGroupSE(ArchRProj = filt_sub_atac_proj, useMatrix = "MotifMatrix", groupBy = "FineNamedClust")
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

# Surfactant production regulators ######
library(biomaRt)

filt_sub_rna_proj <- readRDS(file = "/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Epithelial/filtered_Epithelial.rds")

ensembl=useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')

# Surfactant production related genes
sft_genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'name_1006'), 
      filters = 'go', 
      values = c('GO:0043129', #surfactant homeostasis,
                 "GO:0042599", "GO:0097208", #lamellar body, alveolar lamellar body,
                 "GO:0097233", "GO:0060510" #alveolar lamellar body membrane, type II pneumocyte differentiation
      ),   
      mart = ensembl)

# Add module score
sft_genes <- sft_genes$external_gene_name %>% unique()
sft_genes <- sft_genes[which(sft_genes %in% rownames(filt_sub_rna_proj))]
sft_genes <- list(sft_genes)
filt_sub_rna_proj <- AddModuleScore(filt_sub_rna_proj, features = sft_genes, assay = "RNA", name = "SFTgenes")

# Plot surfactant module expression
p <- FeaturePlot(filt_sub_rna_proj, 
            features = "SFTgenes1", 
            max.cutoff = "q99", 
            min.cutoff = "q1", 
            pt.size = 1, cols = cmaps_BOR$sunrise,
            #coord.fixed = 1,
            order = TRUE #label = TRUE, label.size = 5
            ) + coord_fixed(ratio = 1.5)

pdf(paste0(plotDir, "/UMAP_SFTgene_ModuleScore_sunrise.pdf"))
p + theme_ArchR()
dev.off()

#readr::write_tsv(data.frame(genes = unlist(sft_genes)), paste0(plotDir, "/SFTgenes.tsv"))
sft_genes <- readr::read_tsv(paste0(plotDir, "/SFTgenes.tsv")) %>% unlist() %>% unname()

# Peak module for surfactant module
p2gGR <- getP2G_GR(atac_proj, corrCutoff = 0.45, filtNA = TRUE)
p2gGR_metadata <- mcols(p2gGR)
p2g_sft <- p2gGR[which(p2gGR_metadata$symbol %in% unlist(sft_genes)),]
p2g_sft <- as.list(GRangesList(p2g_sft))

filt_sub_atac_proj <- addModuleScore(filt_sub_atac_proj, 
                                     useMatrix = "PeakMatrix",
                                     name = "SFTgenes_P2G",
                                     features = p2g_sft)

pdf(paste0(plotDir, "/UMAP_SFTgene_Peak2Gene_ModuleScore.pdf"))
plotEmbedding(filt_sub_atac_proj, 
              embedding = "customUMAP", 
              name = "SFTgenes_P2G1", 
              quantCut = c(0.01, 0.99),
              size = 1, rastr = F)
  
dev.off()

# Custom plot of the surfactant peak module
dfToPlot <- data.frame(filt_sub_atac_proj@embeddings$customUMAP$df$`custom#UMAP1`,
                       filt_sub_atac_proj@embeddings$customUMAP$df$`custom#UMAP2`,
                       filt_sub_atac_proj@cellColData$SFTgenes_P2G1)
colnames(dfToPlot) <- c("x", "y", "SFTgenes_P2G")

plotUMAP(dfToPlot, 
         dataType = "quantitative", 
         cmap = cmaps_BOR$solarExtra, point_size = 1,
         colorLims = c(min(dfToPlot$SFTgenes_P2G), max(dfToPlot$SFTgenes_P2G)))

peaks_gr <- getPeakSet(filt_sub_atac_proj)
peaks_gr_near_gene <- peaks_gr[which(peaks_gr$nearestGene %in% unlist(sft_genes)),]
peaks_gr_near_gene <- as.list(GRangesList(peaks_gr_near_gene))

filt_sub_atac_proj <- addModuleScore(filt_sub_atac_proj, 
                                     useMatrix = "PeakMatrix",
                                     name = "peaksNearSFTgenes",
                                     features = peaks_gr_near_gene)

pdf(paste0(plotDir, "/UMAP_SFTgene_PeaksNearGenes_ModuleScore.pdf"))
plotEmbedding(filt_sub_atac_proj, embedding = "customUMAP", name = "peaksNearSFTgenes1", quantCut = c(0.01, 0.99))
dev.off()


# Identify TFs likely to be functional in this compartment
# Correlate motifmatrix and gene expression matrix
alveolar_sub_atac <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Alveolar")

corGE_MM <- correlateMatrices(
  ArchRProj = alveolar_sub_atac,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PCA"
)

corGE_MM <- corGE_MM[order(abs(corGE_MM$cor), decreasing = TRUE), ]
corGE_MM <- corGE_MM[which(!duplicated(gsub("\\-.*","",corGE_MM[,"MotifMatrix_name"]))), ]
regTFs <- corGE_MM$GeneExpressionMatrix_name[which(corGE_MM$cor > 0.1)]

# Filter TFs that have at least a count of 1 in at least one cluster
filt_sub_rna_proj <- subset(filt_sub_rna_proj, idents = c("APr1", "APr2", "APr3", "APr4", "EpiC", "TiP1", "TiP2", "AT2l", "AT1l"))
avg_expr <- AverageExpression(filt_sub_rna_proj, 
                              assays = "RNA", 
                              slot = "counts",
                              group.by = "FineNamedClust", 
                              features = regTFs) %>% as.data.frame()
regTFs <- row.names(avg_expr)[apply(avg_expr, 1, function(row) {any(row > 3)})]

# Remove antisense genes
regTFs <- regTFs[!grepl("-AS1", regTFs)]
saveRDS(regTFs, paste0(plotDir, "regTFs_sft_genes_filtered.rds"))

# Motif enrichments in linked peaks around the sft genes
# Get locations of motifs of interest:
motifPositions <- getPositions(alveolar_sub_atac, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")
motifLocs <- motifGR[motifGR$motifName %in% regTFs] # filter for TF motifs that are possibly regulatory

# Get peak2gene linkages for the target TF
p2gGR <- getP2G_GR(atac_proj, corrCutoff = 0.45, filtNA = TRUE)
p2gGR_metadata <- mcols(p2gGR)
p2g_sft <- p2gGR[which(p2gGR_metadata$symbol %in% unlist(sft_genes)),]

alveolar_sub_atac <- addBgdPeaks(alveolar_sub_atac, force = TRUE)
enrichedMotifs <- customEnrichment(p2g_sft, matches = getMatches(alveolar_sub_atac), bgdPeaks = getBgdPeaks(alveolar_sub_atac))
enrichedMotifs <- filter(enrichedMotifs, feature %in% regTFs) # Keep TFs in the candidate list
enrichedMotifs <- enrichedMotifs[order(-enrichedMotifs$Enrichment), ]

enrichedMotifs %>% head(20)
enrichedMotifs$rank <- seq_len(nrow(enrichedMotifs))

pdf(paste0(plotDir, "/MotifEnrichment_P2G_SFTgenes.pdf"))
ggplot(enrichedMotifs, aes(rank, Enrichment, color = mlog10p)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = enrichedMotifs[rev(seq_len(10)), ], aes(x = rank, y = Enrichment, label = feature), 
    size = 4,
    nudge_x = 1,
    color = "black"
  ) + theme_ArchR() + 
  ylab("Enrichment") + 
  xlab("Rank Sorted TFs") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
dev.off()

# Plot trajectory of gene expression for surfactant genes along AT2 trajectory
trajectory <- c("TiP1", "TiP2", "AT2l")
trajectory_name <- paste(trajectory, collapse = "_")
trajGE <- readRDS(paste0(trajDir, "/TrajGE_", trajectory_name, ".rds"))
rownames(trajGE) <- rownames(trajGE) %>% strsplit(":") %>% sapply('[',2)

# Exclude ribosomal genes
blacklist <- c(
  grep(pattern="^RPS", x=rownames(trajGE), value=TRUE),
  grep(pattern="^RPL", x=rownames(trajGE), value=TRUE)
)
trajGE <- trajGE[rownames(trajGE) %ni% blacklist,]
exclude <- c("NKIRAS1", "SMPD1", "NKIRAS2", "FGF10", "NKX2-1")
trajGE <- trajGE[rownames(trajGE) %ni% exclude,]

varCutoff <- 0.5 # Cutoff for variable genes
varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(assays(trajGE)$smoothMat))
var.traj.genes <- rownames(trajGE[order(varQ, decreasing=TRUE),]) %>% head((1-varCutoff)*nrow(trajGE))
#var.traj.genes <- var.traj.genes[which(var.traj.genes %ni% exclude)]
#include <- c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "LAMP3", "NFIB", "ABCA12",
#             "RAB7A", "LPCAT1", "EPAS1", "ABCA3", "CTSH", "ADGRF5")
#var.traj.genes <- var.traj.genes[which(var.traj.genes %in% include)]
sft.trajGE <- trajGE[var.traj.genes[var.traj.genes %in% unlist(sft_genes)],]

var.cutoff <- 0 # Default is 0.9
lims <- c(-1.5, 1.5) # Defaults are c(-1.5, 1.5)

p1 <- plotTrajectoryHeatmap(sft.trajGE, pal=cmaps_BOR$sunrise, 
                            varCutOff=var.cutoff, limits=lims, labelTop=50)

pdf(paste0(plotDir, "/SFTgene_trajGE_heatmap.pdf"), w = 5, h = 4)
p1
dev.off()

# Browser track plots for linked peaks to sft genes
clustOrder <- fineOrder %>% unlist() %>% unname()
celltypes <- c("TiP1", "TiP2", "AT2l", "AT1l")
genesToPlot <- p2g_sft$symbol %>% unique()

p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = celltypes,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = genesToPlot, 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500
)

plotPDF(plotList = p, 
        name = paste0("SFT_genes_P2G_BrowserPlot.pdf"), 
        ArchRProj = atac_proj, 
        addDOC = FALSE, 
        width = 6, height = 3.5)


####################################
# Find motif positions of the top motifs enriched in surfactant genes
####################################
alveolar_sub_atac <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Alveolar")

# Motif enrichments in linked peaks around the sft genes
# Get locations of motifs of interest:
motifPositions <- getPositions(alveolar_sub_atac, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")
regTFs <- readRDS(paste0(plotDir, "regTFs_sft_genes_filtered.rds"))
motifLocs <- motifGR[motifGR$motifName %in% regTFs] # filter for TF motifs that are possibly regulatory

# Get peak2gene linkages for the target TF
p2gGR <- getP2G_GR(atac_proj, corrCutoff = 0.45, filtNA = TRUE)
p2gGR_metadata <- mcols(p2gGR)

sft_genes <- readr::read_tsv(paste0(plotDir, "/SFTgenes.tsv")) %>% unlist() %>% unname()
p2g_sft <- p2gGR[which(p2gGR_metadata$symbol %in% unlist(sft_genes)),]

# Create bed file for chrombpnet contribution scores

# Get blacklist regions used by chrombpnet
blacklist <- import.bed("/oak/stanford/groups/wjg/skim/projects/LDA/resources/chrombpnet/blacklist.bed.gz")
# Extend blacklist by 1057bp on both sides
blacklist_ext <- resize(blacklist, width = width(blacklist) + (1057 * 2), fix = "center")

filtered_peaks_GR <- p2g_sft[!p2g_sft %over% blacklist_ext]
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
peakDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/Figure04_Epithelial/epi_sft_gene_peaks"
dir.create(peakDir, showWarnings = FALSE, recursive = TRUE)
write.table(bed, file = paste0(peakDir, "/epi_sft_gene_peaks.bed"), quote = F, sep = "\t", row.names = F, col.names = F)

# After running chrombpnet for contribution scores and finemo for hitcalling, identify overlapping NR3C1 motif with hits
motifs_nr3c1 <- motifLocs[motifLocs$motifName == "NR3C1"]
saveRDS(motifLocs, file="motifslocs.rds")
hits <- readr::read_tsv(paste0(wd, "/epi_sft_gene_peaks_hitcalling/AT2l/hits_unique.tsv"))

hits_GR <- GRanges(
  seqnames = hits$chr,
  ranges = IRanges(start = hits$start, end = hits$end),
  strand = hits$strand,
  mcols = hits %>% dplyr::select(-chr, -start, -end, -strand) # Add metadata columns
)

motifs_nr3c1[motifs_nr3c1 %over% hits_GR]
hits_GR[hits_GR %over% motifs_nr3c1]$mcols.motif_name

saveRDS(hits_GR, paste0(wd, "/epi_sft_gene_peaks_hitcalling/nr3c1_over_hits.rds"))

#### Plot TFregs motif deviations ####

# Plot chromVAR deviations for candidate TFs for each vSMC subtype
markerMotifs <- paste0("z:", regTFs)

p <- plotEmbedding(ArchRProj = sub_atac_proj, 
                   groupBy = "FineNamedClust", 
                   colorBy = "MotifMatrix", 
                   name = sort(markerMotifs),
                   embedding = "customUMAP",
                   size = 1
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

pdf(paste0(plotDir, "/surf_regTF_chromVAR_UMAP.pdf"), w = 20, h = 10)
do.call(cowplot::plot_grid, c(list(ncol = 7),p2))
dev.off()


