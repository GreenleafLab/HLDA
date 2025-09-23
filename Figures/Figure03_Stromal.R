#!/usr/bin/env Rscript

########################################
# Figure 3 panels
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
source(paste0(scriptPath, "/sample_metadata.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/Figure03_Stromal"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

# set plot directory
plotDir <- paste0(wd)

# Misc options
addArchRGenome("hg38")
pointSize <- 1

##########################################################################################
# Preparing Data
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
sub_rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Stromal/Stromal.rds")
sub_atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Stromal")

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
# Barplot of proportion of vSMC1 vs. vSMC2 across time points
##########################################################################################
vSMC_rna_proj <- subset(sub_rna_proj, ident = c("vSMC1", "vSMC2"))

ageByVSMC <- fractionXbyY(vSMC_rna_proj$age, vSMC_rna_proj$FineNamedClust, xname = "Gestational Age", yname = "vSMC type", add_total = T)

pdf(paste0(plotDir, "/ageByvSMCStatus_stackedBarPlot.pdf"))
stackedBarPlot(ageByVSMC, barwidth = 1, cmap = FineNamedClustCmap, namedColors = T)
dev.off()

# Also show the distribution of arterial endothelial cells vs venous endothelial cells in the dataset
Idents(rna_proj) <- "FineNamedClust"
vasc_rna_proj <- subset(rna_proj, ident = c("eArtr", "lArtr", "Veno"))
vasc_rna_proj$vasc_ident <- vasc_rna_proj$FineNamedClust
vasc_rna_proj$vasc_ident[which(vasc_rna_proj$FineNamedClust %in% c("eArtr", "lArtr"))] <- gsub("^.", "", vasc_rna_proj$FineNamedClust[which(vasc_rna_proj$FineNamedClust %in% c("eArtr", "lArtr"))])

vasc.cmap <- list("Veno" = "#FEE52C",
                  "Artr" = "#F09FAB"
                  )

ageByVasc <- fractionXbyY(vasc_rna_proj$age, vasc_rna_proj$vasc_ident, xname = "Gestational Age", yname = "Vascular structure", add_total = T)

pdf(paste0(plotDir, "/ageByVascStatus_stackedBarPlot.pdf"))
stackedBarPlot(ageByVasc, barwidth = 1, cmap = vasc.cmap, namedColors = T)
dev.off()

##########################################################################################
# Differential abundance testing for stromal cells
##########################################################################################

library(miloR)

#Set/Create Working Directory to Folder
subPlotDir <- paste0(wd, "/milo")
dir.create(subPlotDir, showWarnings = FALSE, recursive = TRUE)

Idents(rna_proj) <- "compartment"
str_proj <- subset(rna_proj, ident = "Stromal")
str_proj@reductions$umap <- sub_rna_proj@reductions$umap

# Convert to SCE since milo needs it
sce <- as.SingleCellExperiment(str_proj, assay = "SCT")
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
str.order <- c(unlist(unname(fineOrder))[which(unlist(unname(fineOrder)) %in% unique(sub_rna_proj$FineNamedClust))], "Mixed")
da_results$FineNamedClust <- factor(da_results$FineNamedClust, levels = rev(str.order))

# Plot beeswarm
pdf(paste0(subPlotDir, "/miloR_RNA_DA_beeswarm.pdf"), w = 5, h = 3)
plotDAbeeswarm(da_results, group.by = "FineNamedClust") + 
  guides(color=guide_colorbar()) +
  theme_BOR() +
  geom_hline(yintercept = 0)
dev.off()

##########################################################################################
# Differential expression testing between vSMCs
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
pseudo_rna_sub <- subset(pseudo_rna, ident = c("vSMC1", "vSMC2"))
data.use <- GetAssayData(pseudo_rna_sub, "counts", "RNA")

# Add metadata information on sample (batch)
group.info <- data.frame(celltype = as.factor(gsub(".*_", "", unname(colnames(data.use)))),
                         samplename = unname(colnames(data.use)), 
                         sample = gsub("_(.*)", "", unname(colnames(data.use)))
                         )
group.info$celltype <- relevel(group.info$celltype, ref = "vSMC1")

dds1 <- DESeq2::DESeqDataSetFromMatrix(
  countData = data.use,
  colData = group.info,
  design = ~ celltype + sample
)

# Run DESeq2
dds1 <- DESeq(dds1)
#resultsNames(dds1)
resLFC <- lfcShrink(dds1, coef="celltype_vSMC2_vs_vSMC1", type="apeglm")
resOrdered <- resLFC[order(resLFC$pvalue),]

saveRDS(dds1, paste0(DE_dir, "/DDS.rds"))
saveRDS(resOrdered, paste0(DE_dir, "/resOrdered.rds"))
resOrdered <- readRDS(paste0(DE_dir, "/resOrdered.rds"))

res.df <- resOrdered %>% as.data.frame()
res.df$gene <- rownames(res.df)
readr::write_csv(res.df, paste0(DE_dir, "/vSMC_DE_genes.csv"))
readr::write_tsv(res.df, paste0(DE_dir, "/vSMC_DE_genes.tsv"))

# Plot MA Plot of differential GE testing
l2fc_threshold <- 1.0
padj_threshold <- 0.01
geneToLabel <- c("ITGA11", "PLN", "COL1A2", "ELN", "COL3A1", "ADAMTS2", "MEF2C", "COL6A6", "GPC3",
                 "CSMD1", "RCAN2")
res.df$label <- NA
idxToLabel <- which(res.df$gene %in% geneToLabel)
res.df$label[idxToLabel] <- res.df$gene[idxToLabel]

res.df$color <- "grey78"
res.df$color[which(res.df$log2FoldChange > l2fc_threshold & res.df$padj < padj_threshold)] <- "#8A9FD1" #vSMC2
res.df$color[which(res.df$log2FoldChange < -l2fc_threshold & res.df$padj < padj_threshold)] <- "#C06CAB" #vSMC1

res.plot <- res.df[,c("baseMean", "log2FoldChange", "label", "color")]

point_size = 1.5
color_col <- "color"

pdf(paste0(DE_dir, "/MA_plot_vSMC1vsvSMC2_l2fc", l2fc_threshold,"_padj", padj_threshold,".pdf"))
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

# GO term enrichment for diff expr genes
source(paste0(scriptPath, "/GO_wrappers.R"))

fates <- c("vSMC2", "vSMC1")

count.mat <- Seurat::GetAssayData(object=vSMC_rna_proj, slot="counts", assay = "RNA")
minUMIs <- 1
minCells <- 2
all_genes <- rownames(count.mat[rowSums(count.mat > minUMIs) > minCells,])

upGenes <- res.df$gene[which(res.df$log2FoldChange > l2fc_threshold & res.df$padj < padj_threshold)]
downGenes <- res.df$gene[which(res.df$log2FoldChange < -l2fc_threshold & res.df$padj < padj_threshold)]

gene_sets <- list(upGenes, downGenes)
names(gene_sets) <- fates

saveRDS(gene_sets, paste0(DE_dir, "/vSMC_DE_genes.rds"))

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

names(GOresults) <- fates

# Plots of GO term enrichments:
pdf(paste0(DE_dir, sprintf("/vSMC_DE_genes_GO_3termsBPonly_l2fc%s_padj%s.pdf", l2fc_threshold, padj_threshold)), width=10, height=2.5)
for(name in names(GOresults)){
  goRes <- GOresults[[name]]
  if(nrow(goRes)>1){
    print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
                       nterms=3, border_color="black", 
                       barwidth=0.85, title=name, barLimits=c(0, 15)))
  }
}
dev.off()


##########################################################################################
# Differential peak testing between vSMCs
##########################################################################################

markerTest <- getMarkerFeatures(
  ArchRProj = sub_atac_proj, 
  useMatrix = "PeakMatrix",
  groupBy = "FineNamedClust",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "vSMC2",
  bgdGroups = "vSMC1"
)

pma <- plotMarkers(seMarker = markerTest, 
                   name = "vSMC2", 
                   cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", 
                   plotAs = "MA")

pma

motifsvSMC2 <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = sub_atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 1"
)

df <- data.frame(TF = rownames(motifsvSMC2), mlog10Padj = assay(motifsvSMC2)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsvSMC1 <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = sub_atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -1"
)

df <- data.frame(TF = rownames(motifsvSMC1), mlog10Padj = assay(motifsvSMC1)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDown <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf(paste0(DE_dir, "/MA_plot_peaks_vSMC1vsvSMC2_l2fc1_padj0.05.pdf"), w = 5, h = 5)
pma
dev.off()

pdf(paste0(DE_dir, "/vSMC_DE_peaks_MotifEnriched_vSMC2_vSMC1_l2fc1_padj0.05.pdf"), w = 5, h = 5)
ggUp
ggDown
dev.off()

# Differential accessibility on linked peaks
diffGenes <- unlist(gene_sets) %>% unname()

full_p2gGR <- readRDS(file="/oak/stanford/groups/wjg/skim/projects/LDA/Figure02_P2G_Analysis/multilevel_p2gGR_sorted.rds") # NOT merged or correlation filtered

filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()
filt_p2gGR <- filt_p2gGR[!is.na(filt_p2gGR$idxATAC)] # Remove ATAC peak indices that are NA for some reason
# Causes errors with ArchR getPeak2GeneLinks otherwise

# Correlation filter
filt_p2gGR <- filt_p2gGR[which(filt_p2gGR$Correlation > 0.45)]

# Select peaks with linkages to DE genes
diffP2gPeaks <- filt_p2gGR[filt_p2gGR$symbol %in% diffGenes]

subpeaks <- getPeakSet(sub_atac_proj)
filt_subpeaks <- subpeaks[which(subpeaks$nearestGene %in% diffGenes)]

saveArchRProject(sub_atac_proj, paste0(DE_dir, "/stromal_atac_proj_vSMC_DE_genes_peaks_only"), load = FALSE)
vSMC_sub_atac_proj <- loadArchRProject(path = paste0(DE_dir, "/stromal_atac_proj_vSMC_DE_genes_peaks_only"))

vSMC_sub_atac_proj <- addPeakSet(
  ArchRProj = vSMC_sub_atac_proj,
  peakSet = filt_subpeaks,
  force = TRUE
)

getPeakSet(vSMC_sub_atac_proj)

vSMC_sub_atac_proj <- addPeakMatrix(vSMC_sub_atac_proj)

getAvailableMatrices(vSMC_sub_atac_proj)

markerTest.filtGenes <- getMarkerFeatures(
  ArchRProj = vSMC_sub_atac_proj, 
  useMatrix = "PeakMatrix",
  groupBy = "FineNamedClust",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "vSMC2",
  bgdGroups = "vSMC1"
)

pma <- plotMarkers(seMarker = markerTest.filtGenes, 
                   name = "vSMC2", 
                   cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1.0", 
                   plotAs = "MA")

# Load and convert PFMS from FigR paper into PWMs
cisbp_pwms <- readRDS("/oak/stanford/groups/wjg/skim/resources/MotifPWMs/FigR_cisbp_pwms/cisBP_human_pfms_2021.rds")
cisbp_pwms <- TFBSTools::toPWM(cisbp_pwms)
vSMC_sub_atac_proj <- addMotifAnnotations(ArchRProj = vSMC_sub_atac_proj, 
                                          motifPWMs = cisbp_pwms, 
                                          name = "Motif",
                                          force = TRUE)

motifsvSMC2 <- peakAnnoEnrichment(
  seMarker = markerTest.filtGenes,
  ArchRProj = vSMC_sub_atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 1"
)

df <- data.frame(TF = rownames(motifsvSMC2), mlog10Padj = assay(motifsvSMC2)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsvSMC1 <- peakAnnoEnrichment(
  seMarker = markerTest.filtGenes,
  ArchRProj = vSMC_sub_atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -1"
)

df <- data.frame(TF = rownames(motifsvSMC1), mlog10Padj = assay(motifsvSMC1)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDown <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf(paste0(DE_dir, "/MA_plot_peaks_diffGEonly_vSMC1vsvSMC2_l2fc1_padj0.05.pdf"), w = 5, h = 5)
pma
dev.off()

pdf(paste0(DE_dir, "/vSMC_DE_peaks_diffGEonly_MotifEnriched_vSMC2_vSMC1_l2fc1_padj0.05.pdf"), w = 5, h = 5)
ggUp
ggDown
dev.off()

# chromVAR variable deviations
sub_atac_proj_stromal_b <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Stromal_b")
plotVarDev <- getVarDeviations(sub_atac_proj_stromal_b, name = "MotifMatrix", plot = TRUE)

pdf(paste0(DE_dir, "/chromVAR_deviations_stromal_b_only.pdf"))
plotVarDev
dev.off()

# Positive TF regulators
# Identify deviant TF motifs
seGroupMotif <- getGroupSE(ArchRProj = sub_atac_proj_stromal_b, useMatrix = "MotifMatrix", groupBy = "FineNamedClust")
seGroupMotif <- seGroupMotif[,which(colnames(seGroupMotif) %in% c("vSMC1", "vSMC2"))]
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

# Correlate motifmatrix and gene expression matrix
corGE_MM <- correlateMatrices(
  ArchRProj = sub_atac_proj_stromal_b,
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

pdf(paste0(DE_dir, "/PositiveTF_Regulators_stromal_b.pdf"), width = 5, height = 5)
p
dev.off()

DotPlot(sub_rna_proj, features = "AGTR1")

p <- plotBrowserTrack(
  ArchRProj = sub_atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = c("ePeri", "lPeri", "vSMC1", "vSMC2"),
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = c("AGTR1", "MEF2C", "ITGA11", "PLN"), 
  loops = getPeak2GeneLinks(sub_atac_proj),
  tileSize=500
)

plotPDF(plotList = p, 
        name = "/vSMC_interesting_genes_BrowserTracks.pdf", 
        ArchRProj = sub_atac_proj, 
        addDOC = FALSE, 
        width = 6, height = 6)


###### Motif-Gene Correlation by Log2FoldChange Expression between vSMC1 vs. vSMC2 ######

# Correlate motifmatrix and gene expression matrix
sub_atac_proj_stromal_b <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Stromal_b")

corGE_MM <- correlateMatrices(
  ArchRProj = sub_atac_proj_stromal_b,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PCA"
)

corGE_MM <- corGE_MM[order(abs(corGE_MM$cor), decreasing = TRUE), ]
corGE_MM <- corGE_MM[which(!duplicated(gsub("\\-.*","",corGE_MM[,"MotifMatrix_name"]))), ] # remove AS and other non-gene variants of TF

res.df <- readr::read_tsv(paste0(DE_dir, "/vSMC_DE_genes.tsv"))

l2fc_threshold <- 1
padj_threshold <- 0.01
plot.df <- left_join(as.data.frame(corGE_MM), as.data.frame(res.df), by = c("GeneExpressionMatrix_name" = "gene"))
plot.df$label <- "NO"
plot.df$label[which(plot.df$cor > 0.0 & plot.df$padj.x < padj_threshold & abs(plot.df$log2FoldChange) > l2fc_threshold & plot.df$padj.y < padj_threshold)] <- "YES"

pdf(paste0(DE_dir, "/vSMC_Motif_Activity_Expression_Correlation_By_vSMC1vsvSMC2_DiffExp.pdf"), w = 6, h = 6)
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
  ylab("Log2FoldChange vSMC2 vs. vSMC1") +
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(-3, 3), expand = c(0,0)) +
  geom_hline(yintercept = c(-l2fc_threshold, l2fc_threshold), linetype="dashed")
dev.off()

###### Footprints for TFs correlated with each cell type ######
motifPositions <- getPositions(sub_atac_proj_stromal_b)

motifs <- c("NR2F1", "LEF1", "MEF2C", "MEF2B", "TWIST2", "GATA6", "EBF2", "ETS1", "ZNF331", "STAT4", "CREB5", "NFATC2", "PRRX2", "ZNF385D")
#motifs <- plot.df$GeneExpressionMatrix_name[plot.df$label == "YES"]
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

if(is.null(sub_atac_proj_stromal_b@projectMetadata$GroupCoverages$FineNamedClust)){
  sub_atac_proj_stromal_b <- addGroupCoverages(ArchRProj = sub_atac_proj_stromal_b, groupBy = "FineNamedClust")
}

seFoot <- getFootprints(
  ArchRProj = sub_atac_proj_stromal_b, 
  positions = motifPositions[markerMotifs], 
  groupBy = "FineNamedClust", 
  nTop = 10000,
  useGroups = c("vSMC1", "vSMC2", "ePeri", "lPeri")
)

#saveArchRProject(sub_atac_proj_stromal_b, outputDirectory = "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Stromal_b")
#sub_atac_proj_stromal_b <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Stromal_b")

plotFootprints(
  seFoot = seFoot,
  ArchRProj = sub_atac_proj_stromal_b, 
  normMethod = "Subtract",
  plotName = "vSMC-Motifs-Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = sub_atac_proj_stromal_b, 
  normMethod = "Divide",
  plotName = "vSMC-Motifs-Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

##########################################################################################
# Plotting cpdb results
##########################################################################################
library(ktplots)
library(ggplot2)
library(SingleCellExperiment)

# Set/Create Working Directory
results.dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/4_signaling/results_vSMC"

# subset seurat object
subobj <- subset(rna_proj, idents = c("eArtr", "lArtr", "Veno", "vSMC1", "vSMC2"))
subobj$FineNamedClust[which(subobj$FineNamedClust %in% c("eArtr", "lArtr"))] <- "Artr"

Idents(subobj) <- "FineNamedClust"
levels(subobj) <- c("Artr", "Veno", "vSMC1", "vSMC2")

scdata <- as.SingleCellExperiment(subobj, assay = "RNA")

pvals <- read.delim(paste0(results.dir, "/degs_analysis_relevant_interactions_02_27_2024_151501.txt"), check.names = FALSE)
means <- read.delim(paste0(results.dir, "/degs_analysis_means_02_27_2024_151501.txt"), check.names = FALSE)
interaction <- read.delim(paste0(results.dir, "/degs_analysis_interaction_scores_02_27_2024_151501.txt"), check.names = FALSE)
deconvoluted <- read.delim(paste0(results.dir, "/degs_analysis_deconvoluted_02_27_2024_151501.txt"), check.names = FALSE)

plot_cpdb_heatmap(pvals = pvals, 
                  degs_analysis = TRUE,
                  title = "Sum of significant interactions")

# All significant interactiosn with alpha based on interactio nscore (specificity)
p <- plot_cpdb(cell_type1 = 'Artr|Veno', 
          cell_type2 = 'vSMC1|vSMC2', 
          scdata = subobj,
          celltype_key = "FineNamedClust",
          means = means, 
          pvals = pvals,
          #genes = c("DLL4", "RSPO3", "NTF3"),
          highlight_size=1,
          degs_analysis=TRUE,
          standard_scale=TRUE,
          interaction_scores = interaction,
          scale_alpha_by_interaction_scores=TRUE, 
          
          #min_interaction_score = 20
) + small_guide() + small_legend(keysize=1) 

pdf(paste0(plotDir, "/cpdb_plot_AllInteractions.pdf"), w = 20, h = 200)
p + facet_wrap(~classification, ncol = 5)
dev.off()

Idents(rna_proj) <- "FineNamedClust"
levels(rna_proj) <- fineOrder %>% unname() %>% unlist()
DotPlot(rna_proj, features = c("WNT2B", "SFRP1", "FRZB", "FZD2", "LRP6"))
DotPlot(rna_proj, features = c("DLL1", "DLL4", "JAG1" ,"JAG2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"))
DotPlot(rna_proj, features = c("EFNB2", "EFNB1", "EFNA1", "EPHA3", "EPHA4", "EPHA7", "EPHA6"))

pdf(paste0(plotDir, "/cpdb_plot__SignificantInteractions.pdf"), w = 6, h = 6)
plot_cpdb(cell_type1 = 'Artr|Veno', 
          cell_type2 = 'vSMC1|vSMC2', 
          scdata = subobj,
          celltype_key = "FineNamedClust",
          means = means, 
          pvals = pvals,
          genes = c("DLL4", "RSPO3", "WNT2B", "VEGFA", "EFNB2"),
          highlight_size=1,
          degs_analysis=TRUE,
          standard_scale=TRUE,
          interaction_scores = interaction,
          scale_alpha_by_interaction_scores=TRUE,
          col_option=cmaps_BOR$sunrise) + 
  small_guide() + 
  small_legend(keysize=1)
dev.off()

ligands <- c("VEGFA", "SFRP1", "FZD2", "FRZB", "EPHB6", "RSPO3", "DLL4")
subclustOrder <- c("Artr", "Veno", "vSMC1", "vSMC2")
levels(subobj) <- subclustOrder
pdf(paste0(plotDir, "/cpdb_dotplot_ligands.pdf"), w = 5, h = 5)
DotPlot(subobj, features = ligands) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  scale_color_gradientn(colors = cmaps_BOR$whitePurple) +
  geom_point(aes(size=pct.exp), shape = 21, colour="white", stroke=0.5) +
  theme_BOR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()

receptors <- c("NRP2", "KDR", "FLT1", "WNT2B", "EFNB2", "LGR4", "NOTCH3", "NOTCH2")
#subclustOrder <- c("vSMC1", "vSMC2", "Artr", "Veno")  
#levels(subobj) <- subclustOrder
pdf(paste0(plotDir, "/cpdb_dotplot_receptors.pdf"), w = 5, h = 5)
DotPlot(subobj, features = receptors) +
  scale_x_discrete(position = "bottom") +
  coord_flip() +
  scale_color_gradientn(colors = cmaps_BOR$whitePurple) +
  geom_point(aes(size=pct.exp), shape = 21, colour="white", stroke=0.5) +
  theme_BOR() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()

# Norm count dotplots ####
aspectRatio <- 2

# Dot plot of ligands
count_mat <- GetAssayData(object=subobj, slot="counts")
avgPctMat <- avgAndPctExpressed(count_mat, subobj$FineNamedClust, feature_normalize=TRUE, min_pct=0)
# Subset to genes we care about:
avgPctMat <- avgPctMat[avgPctMat$feature %in% ligands,]
avgPctMat <- avgPctMat[avgPctMat$grp %in% subclustOrder,]
# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0
# Scale percent expression to max per gene
avgPctMat <- avgPctMat %>% group_by(feature) %>% mutate(scaledPctExpr = pctExpr / max(pctExpr)) %>%
  ungroup() %>% as.data.frame()
# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")

pdf(paste0(plotDir, "/cpdb_dotplot_scaledPctExpr_ligands.pdf"), width=5, height=5)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="scaledPctExpr", 
        xorder=subclustOrder, yorder=ligands, cmap=BuenColors::jdb_palette("brewer_heat", n = 7, type ="discrete"), aspectRatio=aspectRatio) +
  scale_y_discrete(position = "right")
dev.off()

# Dot plot of receptors
avgPctMat <- avgAndPctExpressed(count_mat, subobj$FineNamedClust, feature_normalize=TRUE, min_pct=0)
# Subset to genes we care about:
avgPctMat <- avgPctMat[avgPctMat$feature %in% receptors,]
avgPctMat <- avgPctMat[avgPctMat$grp %in% subclustOrder,]
# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0
# Scale percent expression to max per gene
avgPctMat <- avgPctMat %>% group_by(feature) %>% mutate(scaledPctExpr = pctExpr / max(pctExpr)) %>%
  ungroup() %>% as.data.frame()
# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col="feature", col_col="grp", val_col="avgExpr")

pdf(paste0(plotDir, "/cpdb_dotplot_scaledPctExpr_receptors.pdf"), width=5, height=5)
dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="scaledPctExpr", 
        xorder=subclustOrder, yorder=receptors, cmap=BuenColors::jdb_palette("brewer_blue", n = 7, type ="discrete"), aspectRatio=aspectRatio)
dev.off()


##########################################################################################
# Marker TFs for each of the stromal cell types
##########################################################################################

str_markers <- FindAllMarkers(sub_rna_proj, assay = "RNA", logfc.threshold = 0.5, min.diff.pct = 0.3, only.pos = T)
saveRDS(str_markers, paste0(wd, "/Stromal_RNA_markergenes_logfc0.5_mindiffpct0.3.rds"))
#str_markers <- readRDS(paste0(wd, "/Stromal_RNA_markergenes_logfc0.5_mindiffpct0.3.rds"))

corGE_MM <- correlateMatrices(
  ArchRProj = sub_atac_proj,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PCA"
)

corGE_MM <- corGE_MM[order(abs(corGE_MM$cor), decreasing = TRUE), ]
corGE_MM <- corGE_MM[which(!duplicated(gsub("\\-.*","",corGE_MM[,"MotifMatrix_name"]))), ] # remove AS and other non-gene variants of TF

log10pCut <- 10
l2fc_threshold <- 0.5
padj_threshold <- 0.01
plot.df <- inner_join(as.data.frame(corGE_MM), as.data.frame(str_markers), by = c("GeneExpressionMatrix_name" = "gene"))
plot.df$label <- "NO"
plot.df$label[which(plot.df$cor > 0 & plot.df$padj < padj_threshold & abs(plot.df$avg_log2FC) > l2fc_threshold & plot.df$p_val_adj < padj_threshold)] <- "YES"

plot.df.filt <- plot.df %>% dplyr::filter(label == "YES")


compartment <- "Stromal"
enrichMotifs <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/Figure01_Data_Summary/enrichMotifs-Stromal.rds")
ctToKeep <- fineOrder[[compartment]]
ctToKeep <- ctToKeep[ctToKeep %ni% c("Meso", "Schw")]
plot_mat <- plotEnrichHeatmap(enrichMotifs[unique(plot.df.filt$GeneExpressionMatrix_name),ctToKeep], n=10, transpose=FALSE, 
                              cutOff=log10pCut, returnMatrix=TRUE)

plot_mat <- prettyOrderMat(plot_mat[,ctToKeep], clusterCols=FALSE, cutOff=1)$mat

pdf(paste0(plotDir, "/MarkerPeak-MotifEnriched-corTF-DiffGE-Heatmap-", compartment,".pdf"), width=10, height=14)
fontsize <- 8
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=ctToKeep,col=list(atac_cluster=FineNamedClustCmap), 
                        show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  top_annotation = ta,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = fontsize),
  column_names_gp = gpar(fontsize = fontsize),
  width = ncol(plot_mat)*unit(0.3, "cm"),
  height = nrow(plot_mat)*unit(0.3, "cm"),
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]",
  border_gp = gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()

# Plot chromVAR deviations for candidate TFs for each vSMC subtype
motifNames <- c("NR2F1", "TWIST2", "GATA6", "ETS1", "STAT4", "CREB5", "NFATC2", "PRRX2", "LEF1" ,"MEF2C", "MEF2B", "EBF2")
markerMotifs <- paste0("z:", motifNames)

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
pdf(paste0(plotDir, "/corTF-DiffGE-TFs_chromVAR_UMAP.pdf"), w = 8, h = 6)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
dev.off()

pdf(paste0(plotDir, "/corTF-DiffGE-TFs_GE_UMAP.pdf"), w = 16, h = 10)
FeaturePlot(sub_rna_proj, features = motifNames, pt.size = 0.5, cols = cmaps_BOR$blueYellow)
dev.off()

