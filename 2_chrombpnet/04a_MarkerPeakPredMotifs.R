#!/usr/bin/env Rscript

########################################
# Predict motifs in marker peaks
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
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/markerpeak_motifs_pred"

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
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds")

# Load ArchR Project
atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda")
rna_obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

##########################################################################################
# Export bigwigs for all of the clusters
##########################################################################################
#coverageFiles <- getGroupBW(ArchRProj = atac_proj,
#                            groupBy = "FineNamedClust")


##########################################################################################
# Prepare BED files of marker peaks
##########################################################################################

markerPeaks <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_markerPeaks_FineNamedClust.rds")
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T) # return as a list of GRanges

# Get blacklist regions used by chrombpnet
blacklist <- import.bed("/oak/stanford/groups/wjg/skim/projects/LDA/resources/chrombpnet/blacklist.bed.gz")
# Extend blacklist by 1057bp on both sides
blacklist_ext <- resize(blacklist, width = width(blacklist) + (1057 * 2), fix = "center")

# Create dir for storing bed files
peaks.dir <- paste0(wd, "/peaks")
dir.create(peaks.dir, showWarnings = TRUE, recursive = TRUE)

# Create bed file for each cluster
clusters <- names(markerList)

for (cluster in clusters) {
  # Grab marker peaks for each cluster and remove blacklisted regions
  granges <- markerList[[cluster]]
  granges <- granges[!granges %over% blacklist_ext]

  # skip this cluster if there are no markerpeaks (i.e. Schwann cells don't pass this filter)
  if (length(granges) == 0) {
    next
  }

  # Give name to peaks
  granges$peakName <- (granges %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})

  # add seqinfo
  seqinfo(granges) <- seqinfo(genome)[seqlevels(granges)]

  # Create BED file for each cluster
  bed <- data.frame(
    seqnames = seqnames(granges),
    starts = as.integer(start(granges) - 251), #BED files are 0-indexed but granges are not
    ends = as.integer(end(granges) + 250), # Convert to integer to prevent bed files being written in sci notation
    names = granges$peakName,
    scores = granges$Log2FC, #This score is NOT used by chrombpnet but adding Log2FC for markerpeaks to hold the place
    strands = ".",
    col_6 = ".",
    col_7 = ".",
    col_8 = ".",
    summit = 500 #chrombpnet uses the 9th column as the summit position to calculate gc content. 500 since ArchR creates fixed width peaks that are 501bp long
    # and chrombpnet takese in 1000bp peak regions based on Lakshman paper
  )
  write.table(bed, file = paste0(peaks.dir, "/", cluster, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
}
























##########################################################################################
# PROX1 enahncer
##########################################################################################

# Enhancer for PROX1 described in https://www.nature.com/articles/s41586-022-05650-9/figures/1

prox1_enh <- GRanges(seqnames = "chr1",
        ranges = IRanges(213977514:213977842))

# Check peaks at PROX1 promoter
p <- plotBrowserTrack(
  ArchRProj = atac_proj,
  groupBy = "CustomNamedClust",
  useGroups = c("egCap", "gCap", "cycgCap", "eAero" ,"Aero", "Lymp", "ArtrTip", "Artr1", "Artr2", "Veno", "PNEC1", "PNEC2", "PNEC3"),
  geneSymbol = "PROX1",
  upstream = 35000,
  downstream = 50000,
  loops = getCoAccessibility(atac_proj),
  highlight = prox1_enh
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

# Get PROX1 enhancer peak 
peaks <- getClusterPeaks(atac_proj, clusterNames = "Lymp", peakGR = peaks_gr, groupBy = "CustomNamedClust")
peaks <- peaks[!peaks %over% blacklist_ext]
peaks <- peaks[peaks %over% prox1_enh]

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


matched_motifs_positions <- motif_positions[lengths(motif_positions) != 0]

hits1 <- findOverlaps(motif_gr1, matched_motifs_positions)
matched_motifs_positions1 <- matched_motifs_positions[hits1@to]

hits2 <- findOverlaps(motif_gr2, matched_motifs_positions)
matched_motifs_positions2 <- matched_motifs_positions[hits2@to]

##################################################
# Violin plots of RNA expression for select genes
##################################################
geneToPlot <- c("PROX1", "NFATC2")
clustToPlot <- c("egCap", "gCap", "cycgCap", "eAero" ,"Aero", "Lymp", "ArtrTip", "Artr1", "Artr2", "Veno", "PNEC1", "PNEC2", "PNEC3")

subobj <- rna_obj[,which(rna_obj$CustomNamedClust %in% clustToPlot)]
Idents(subobj) <- "CustomNamedClust"
subobj@active.ident <- factor(x = subobj@active.ident, levels = clustToPlot)

pdf(paste0(plotDir, "/PROX1_Expression_VlnPlot.pdf"), width=10, height=3)
VlnPlot(subobj, features = "PROX1", pt.size = 0) + 
  scale_fill_manual(values = CustomNamedClustCmap[clustToPlot]) +
  xlab("") + 
  theme(legend.position = "none") +
  scale_y_continuous(position = "right")
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
