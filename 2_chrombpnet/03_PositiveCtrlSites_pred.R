#!/usr/bin/env Rscript

########################################
# Predict positive ctrl sites
########################################

# Load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
})


# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/pos_ctrl_pred"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = TRUE, recursive = TRUE)
setwd(wd)

plotDir <- paste0(wd, "/plots")
dir.create(plotDir, showWarnings = TRUE, recursive = TRUE)

# Misc options
addArchRGenome("hg38")

# Load genome
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

# Block scientific notations
options(scipen = 999)

# Load color maps
compartmentCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/lungClusterColors.rds") %>% unlist()
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds") %>% unlist()
sample_cmap <- readRDS(paste0("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds"))
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")

# Load data objects
atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda")
rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# Load subprojects
sub_rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final/Endothelial/Endothelial.rds")
sub_atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/Endothelial")

Idents(sub_rna_proj) <- "FineNamedClust"

##########################################################################################
# PROX1 enahncer
##########################################################################################

# Enhancer for PROX1 described in https://www.nature.com/articles/s41586-022-05650-9/figures/1

prox1_enh <- GRanges(seqnames = "chr1",
        ranges = IRanges(213977514:213977842))

# Check peaks at PROX1 promoter
p <- plotBrowserTrack(
  ArchRProj = atac_proj,
  groupBy = "FineNamedClust",
  useGroups = c("egCap", "lgCap", "eAero" ,"lAero", "eArtr", "lArtr", "Veno", "Lymp", "PNEC1", "PNEC2", "PNEC3"),
  pal = FineNamedClustCmap,
  geneSymbol = "PROX1",
  upstream = 35000,
  downstream = 50000,
  loops = getPeak2GeneLinks(atac_proj),
  highlight = prox1_enh,
  tileSize = 500
)

pdf(paste0(plotDir, "/PROX1_Enhancer_BrowserTrack.pdf"), w = 5, h = 5)
grid.draw(p$PROX1)
dev.off()

# Save peak that overlaps with this PROX1 promoter

# Get all peaks from ArchR project
peaks_gr <- getPeakSet(atac_proj)
peaks_gr$peakName <- (peaks_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
seqinfo(peaks_gr) <- seqinfo(genome)[seqlevels(peaks_gr)]
# Get blacklist regions used by chrombpnet
blacklist <- import.bed("/oak/stanford/groups/wjg/skim/projects/LDA/resources/chrombpnet/blacklist.bed.gz")
# Extend blacklist by 1057bp on both sides
blacklist_ext <- resize(blacklist, width = width(blacklist) + (1057 * 2), fix = "center")

# Pull all linked peaks for PROX1
p2g_GR <- getPeak2GeneLinks(atac_proj, returnLoops = F)
idxATAC <- dplyr::filter(as.data.frame(p2g_GR), 
                         idxRNA == which(metadata(p2g_GR)$geneSet$name == "PROX1")) %>% 
  dplyr::select(idxATAC) %>% 
  as.list() %>% unlist() %>% unname()
p2g_peaks <- p2g_GR@metadata$peakSet[idxATAC]

# Get PROX1 candidate enhancer peaks
peaks <- getClusterPeaks(atac_proj, clusterNames = "Lymp", peakGR = peaks_gr, groupBy = "FineNamedClust")
peaks <- peaks[!peaks %over% blacklist_ext]
peaks <- peaks[peaks %over% c(prox1_enh, p2g_peaks)]

# Create BED file
bed <- data.frame(
  seqnames = seqnames(peaks),
  starts = as.integer(start(peaks) - 251), #BED files are 0-indexed but granges are not
  ends = as.integer(end(peaks) + 250), # Conver to integer to prevent bed files being written in sci notation
  names = peaks$peakName,
  scores = peaks$score,
  strands = ".",
  col_6 = ".",
  col_7 = ".",
  col_8 = ".",
  summit = 500 #chrombpnet uses the 9th column as the summit position to calculate gc content. 500 since ArchR creates fixed width peaks that are 501bp long
  # and chrombpnet takese in 1000bp peak regions based on Laksshman paper
)
write.table(bed, file = paste0(wd, "/PROX1_enh_peak.bed"), quote = F, sep = "\t", row.names = F, col.names = F)


#################################
# Motif matching from contribution scores
#################################

# Match motifs for this region
library(motifmatchr)
library(GenomicRanges)

# load some example motifs
motifs <- readRDS("/oak/stanford/groups/wjg/skim/resources/MotifPWMs/FigR_cisbp_pwms/cisBP_human_pfms_2021.rds")

# Get motif matches in peak
motif_ix <- matchMotifs(motifs, peaks, genome = "hg38") 
matches <- motifMatches(motif_ix) # Extract matches matrix from SummarizedExperiment result
matches_df <- as.matrix(matches) %>% as.data.frame()
matches_df$celltype <- rownames(matches_df)
matches_df <- matches_df %>% pivot_longer(cols = !celltype, names_to = "TF", values_to = "motifMatched")

matched_motifs <- matches_df$TF[which(matches_df$motifMatched == T)]

# Get motif positions within peaks for example motifs in peaks 
motif_positions <- matchMotifs(motifs, peaks, genome = "hg38",
                        out = "positions") 

motif_gr1 <- GRanges(seqnames = "chr1",
                    ranges = IRanges(213977661:213977667))

motif_gr2 <- GRanges(seqnames = "chr1",
                     ranges = IRanges(213977668:213977677))

motif_gr3 <- GRanges(seqnames = "chr1",
                     ranges = IRanges(213977100:213977114)) # Unique to PNECs

matched_motifs_positions <- motif_positions[lengths(motif_positions) != 0]

hits1 <- findOverlaps(motif_gr1, matched_motifs_positions)
matched_motifs_positions1 <- matched_motifs_positions[hits1@to]

hits2 <- findOverlaps(motif_gr2, matched_motifs_positions)
matched_motifs_positions2 <- matched_motifs_positions[hits2@to]

hits3 <- findOverlaps(motif_gr3, matched_motifs_positions)
matched_motifs_positions3 <- matched_motifs_positions[hits3@to]

##################################################
# Violin plots of RNA expression for select genes
##################################################
geneToPlot <- c("PROX1", "NFATC2", "NFATC1", "HNF1A", "HNF1B")
clustToPlot <- c("egCap", "lgCap", "eAero" ,"lAero", "eArtr", "lArtr", "Veno", "Lymp", "PNEC1", "PNEC2", "PNEC3")

subobj <- rna_proj[,which(rna_proj$FineNamedClust %in% clustToPlot)]
Idents(subobj) <- "FineNamedClust"
subobj@active.ident <- factor(x = subobj@active.ident, levels = clustToPlot)

pdf(paste0(plotDir, "/PosCtrlSites_GeneExpression_VlnPlot.pdf"), width=10, height=3)
for (gene in geneToPlot){
  p <- VlnPlot(subobj, features = gene, pt.size = 0) + 
    scale_fill_manual(values = FineNamedClustCmap[clustToPlot]) +
    xlab("") + 
    theme(legend.position = "none") +
    scale_y_continuous(position = "right")
  print(p)
}
dev.off()

label_genes <- c("PROX1", "NFATC2", "NFATC1", "HNF1A", "HNF1B")
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

clustsToPlot <- c("egCap", "lgCap", "eAero" ,"lAero", "eArtr", "lArtr", "Veno", "Lymp", "PNEC1", "PNEC2", "PNEC3")

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

pdf(paste0(plotDir, "/PosCtrlSites_GeneExpression_VlnPlot_ExpressionMat.pdf"), width=10, height=3)
pList
dev.off()






# Plot the expression of the most likely TF binding sites
DotPlot(rna_obj, features = c("GATA2","NFATC2", "NFATC3", "NFATC4", "FOXC2"))

sftpb_cts <- c("APr1", "APr2", "eAlvPr", "AlvPr", "AT1l", "AT2l", "Cili")

genes <- c("SFTPA1","SFTPA2", "SFTPB", "SFTPC", "SFTPD", "ABCA3", "LAMP3")

p <- plotBrowserTrack(
  ArchRProj = atac_proj,
  groupBy = "CustomNamedClust",
  useGroups = sftpb_cts,
  geneSymbol = genes,
  upstream = 100000,
  downstream = 100000,
  loops = getCoAccessibility(atac_proj)
)

pdf(paste0(plotDir, "/surfactant_genes_accessibility.pdf"), w = 5, h = 5)
for (gene in genes) {
  plot <- p[gene]
  grid.draw(plot[[1]])
  grid.newpage()
}
dev.off()

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)

grid.draw(p$LAMP3)


#
