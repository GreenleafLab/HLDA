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
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
httr::set_config(httr::config(ssl_verifypeer = FALSE))
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/Figure02_P2G_Analysis"

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

# P2G definition cutoffs
corrCutoff <- 0.45       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.5  # Default in getCoAccessibility is 0.5

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

##########################################################################################
# Prepare full-project peak to gene linkages, loops, and coaccessibility (full and subproject links)
##########################################################################################
# Prepare lists to store peaks, p2g links, loops, coaccessibility
plot_loop_list <- list()
plot_loop_list[["lda"]] <- getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff, resolution = 1)[[1]]

coaccessibility_list <- list()
coAccPeaks <- getCoAccessibility(atac_proj, corCutOff=corrCutoff, returnLoops=TRUE)[[1]]
coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
coAccPeaks$source <- "lda"
coaccessibility_list[["lda"]] <- coAccPeaks

peak2gene_list <- list()
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=TRUE)
p2gGR$source <- "lda"
peak2gene_list[["lda"]] <- p2gGR

# # Get distribution of peaks to gene linkages and identify 'highly-regulated' genes
# p2gFreqs <- getFreqs(p2gGR$symbol)
# x <- 1:length(p2gFreqs)
# rank_df <- data.frame(npeaks=p2gFreqs, rank=x)

# p <- plotBrowserTrack(
#   ArchRProj = atac_proj,
#   groupBy = "FineNamedClust",
#   geneSymbol = c("NOTCH1", "RUNX3", "TBX4", "MYH11", "CD93", "FOXF1"),
#   upstream = 50000,
#   downstream = 50000,
#   loops = getPeak2GeneLinks(atac_proj)
# )
# 
# grid::grid.draw(p$FOXF1)

# Retrieve information from subclustered objects
subclustered_projects <- c("Epithelial", "Endothelial", "Stromal") # Immune removed due to spurious P2G links

for(subgroup in subclustered_projects){
  message(sprintf("Reading in subcluster %s", subgroup))
  # Read in subclustered project
  sub_dir <- sprintf("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/%s", subgroup)
  sub_proj <- loadArchRProject(sub_dir)
  
  # Get sub-project p2g links
  subP2G <- getP2G_GR(sub_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
  subP2G$source <- subgroup
  peak2gene_list[[subgroup]] <- subP2G
  
  # Get sub-project loops
  plot_loop_list[[subgroup]] <- getPeak2GeneLinks(sub_proj, corCutOff=corrCutoff, resolution = 1)[[1]]
  
  # Get coaccessibility
  coAccPeaks <- getCoAccessibility(sub_proj, corCutOff=coAccCorrCutoff, returnLoops=TRUE)[[1]]
  coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
  coAccPeaks$source <- subgroup
  coaccessibility_list[[subgroup]] <- coAccPeaks
}

full_p2gGR <- as(peak2gene_list, "GRangesList") %>% unlist()
full_coaccessibility <- as(coaccessibility_list, "GRangesList") %>% unlist()

# Fix idxATAC to match the full peak set
idxATAC <- peak2gene_list[["lda"]]$idxATAC
names(idxATAC) <- peak2gene_list[["lda"]]$peakName
full_p2gGR$idxATAC <- idxATAC[full_p2gGR$peakName]

# Save lists of p2g objects, etc.
saveRDS(full_p2gGR, file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
saveRDS(full_coaccessibility, file=paste0(wd, "/multilevel_coaccessibility.rds"))
saveRDS(plot_loop_list, file=paste0(wd, "/multilevel_plot_loops.rds"))

# Load lists of p2g objects, etc.
full_p2gGR <- readRDS(file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
full_coaccessibility <- readRDS(file=paste0(wd, "/multilevel_coaccessibility.rds"))
plot_loop_list <- readRDS(file=paste0(wd, "/multilevel_plot_loops.rds"))

##########################################################################################
# Upset plot of number of peak to gene links identified per group
##########################################################################################

library(UpSetR)

groups <- c("lda", "Epithelial", "Stromal", "Endothelial")
sub_p2gGR <- full_p2gGR[!is.na(full_p2gGR$Correlation)]

upset_list <- lapply(groups, function(g){
  gr <- sub_p2gGR[sub_p2gGR$source == g & sub_p2gGR$Correlation > corrCutoff & sub_p2gGR$VarQATAC > varCutoffATAC & sub_p2gGR$VarQRNA > varCutoffRNA]
  paste0(gr$peakName, "-", gr$symbol)
})
names(upset_list) <- groups


plotUpset <- function(plist, main.bar.color="red", keep.order=FALSE){
  # Function to plot upset plots
  upset(
    fromList(plist), 
    sets=names(plist), 
    order.by="freq",
    #empty.intersections="off",
    point.size=3, line.size=1.5, matrix.dot.alpha=0.5,
    main.bar.color=main.bar.color,
    keep.order=keep.order
  )
}

pdf(paste0(plotDir, "/upset_plot_linked_peaks.pdf"), width=10, height=5)
plotUpset(upset_list, main.bar.color="royalblue1", keep.order=FALSE)
dev.off()

##########################################################################################
# Filter redundant peak to gene links
##########################################################################################

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]

# Save ordered full_p2gGR
saveRDS(full_p2gGR, file=paste0(wd, "/multilevel_p2gGR_sorted.rds")) # NOT merged or correlation filtered
#full_p2gGR <- readRDS(file=paste0(wd, "/multilevel_p2gGR_sorted.rds"))

filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()
filt_p2gGR <- filt_p2gGR[!is.na(filt_p2gGR$idxATAC)] # Remove ATAC peak indices that are NA for some reason
# Causes errors with ArchR getPeak2GeneLinks otherwise

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF


########################################################################################## 
# Compare between linked peaks and "unlinked" peaks
##########################################################################################

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)

all_linked_peaks <- p2gGR$peakName %>% unique()
unlinked_peaks <- allPeaksGR[allPeaksGR$peakName %ni% all_linked_peaks]$peakName

df <- data.frame(
  peakName=c(all_linked_peaks, unlinked_peaks),
  label=c(rep("linked", length(all_linked_peaks)), rep("unlinked", length(unlinked_peaks)))
)

df$GC <- allPeaksGR[df$peakName]$GC
df$log10GeneDist <- log10(allPeaksGR[df$peakName]$distToGeneStart)

# Violin plot of GC content between types
wtest <- wilcox.test(df$GC[df$label == "linked"], df$GC[df$label == "unlinked"])
lmedian <- median(df$GC[df$label == "linked"])
ulmedian <- median(df$GC[df$label == "unlinked"])
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)
cmap <- getColorMap(cmaps_BOR$sambaNight, n=2, type="quantitative")
p <- (
  ggplot(df, aes(x=label, y=GC, fill=label))
  + geom_violin(aes(fill=label), adjust = 1.0, scale='width', position=dodge)
  + stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..),
                 width=0.75, position=dodge, show.legend = FALSE)
  + scale_color_manual(values=cmap)
  + scale_fill_manual(values=cmap)
  + guides(fill=guide_legend(title=""), 
           colour=guide_legend(override.aes = list(size=5)))
  + ggtitle(sprintf("Wilcoxon test pval: %s\n linked median: %s, unlinked median: %s", wtest$p.value, lmedian, ulmedian))
  + xlab("")
  + ylab("GC content")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          aspect.ratio = 1.67,
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1)) 
)
pdf(paste0(plotDir, "/linked_vs_unlinked_GC_content.pdf"), width=6, height=6)
p
dev.off()

# ECDF plot of distance to nearest TSS
p <- (
  ggplot(df, aes(x=log10GeneDist, color=label))
  + stat_ecdf(size=2)
  + geom_vline(xintercept=log10(250000), linetype="dashed") # Dashed line indicating longest possible p2g link
  + scale_color_manual(values=cmap)
  + xlab("Distance to TSS")
  + ylab("Fraction of peaks")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 1)) 
)
pdf(paste0(plotDir, "/linked_vs_unlinked_distance_to_TSS.pdf"), width=6, height=5)
p
dev.off()

##########################################################################################
# Conservation of linked peaks vs unlinked peaks
##########################################################################################
library(phastCons100way.UCSC.hg38)
phast <- phastCons100way.UCSC.hg38

filt_full_p2gGR <- full_p2gGR[full_p2gGR$Correlation > corrCutoff & 
                                full_p2gGR$VarQATAC > varCutoffATAC & 
                                full_p2gGR$VarQRNA > varCutoffRNA]

filt_full_p2gGR$p2gName <- paste0(filt_full_p2gGR$peakName, "_", filt_full_p2gGR$symbol)

# For each source, get the phastCons100 way mean conservation in linked peaks
allPeaksCons <- gscores(phast, allPeaksGR, summaryFun=mean) # Mean is the default summaryFun

# Get conservation of linked and unlinked peaks from each dataset
p2g_groups <- unique(filt_p2gGR$source)

cons_df <- lapply(p2g_groups, function(g){
  sub_p2g_names <- filt_full_p2gGR$peakName[filt_full_p2gGR$source == g] %>% unique()
  lcons <- allPeaksCons$default[allPeaksCons$peakName %in% sub_p2g_names]
  data.frame(
    group=rep(g, times=length(lcons)), 
    conservation=lcons
  )
}) %>% do.call(rbind,.)

# Add conservation of peaks with no link
ulcons <- allPeaksCons$default[allPeaksCons$peakName %ni% filt_full_p2gGR$peakName]
cons_df <- rbind(cons_df, data.frame(
  group=rep("unlinked", times=length(ulcons)),
  conservation=ulcons
))

cons_pvals <- lapply(p2g_groups, function(g){
  glcons <- cons_df[cons_df$group == g, "conservation"]
  ulcons <- cons_df[cons_df$group == "unlinked", "conservation"]
  c(linked_mean=mean(glcons, na.rm=TRUE), unlinked_mean=mean(ulcons, na.rm=TRUE), pval=wilcox.test(glcons, ulcons)$p.value)
}) %>% do.call(rbind,.) %>% as.data.frame()
rownames(cons_pvals) <- p2g_groups

cons_cmap <- rep("royalblue1", times=length(p2g_groups))
names(cons_cmap) <- p2g_groups
cons_cmap <- c(cons_cmap, unlinked="grey")

# Order by decreasing mean conservation
cons_order <- rownames(cons_pvals[order(cons_pvals$linked_mean, decreasing=TRUE),])
cons_df$group <- factor(cons_df$group, levels=c(cons_order, "unlinked"), ordered=TRUE)

p <- (
  ggplot(cons_df, aes(x=group, y=conservation, fill=group), color="black")
  + geom_boxplot(alpha=1.0)
  + scale_y_continuous(limits=c(0.0,1.05), expand=c(0,0))
  + scale_color_manual(values=cons_cmap)
  + scale_fill_manual(values=cons_cmap)
  + xlab("")
  + ylab("phastCons.100")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          legend.position="none", # Remove legend
          axis.text.x = element_text(angle=90, hjust=1)) 
)

pdf(paste0(plotDir, "/linked_vs_unlinked_peak_conservation.pdf"), width=4, height=5)
p
dev.off()

##########################################################################################
# Identify 'highly regulated' genes and 'highly-regulating' peaks
##########################################################################################

# Get all expressed genes:
count.mat <- Seurat::GetAssayData(object=rna_proj, slot="counts")
minUMIs <- 2
minCells <- 3
valid.genes <- rownames(count.mat[rowSums(count.mat > minUMIs) > minCells,])

# Get distribution of peaks to gene linkages and identify 'highly-regulated' genes
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)
p2gFreqs <- getFreqs(p2gGR$symbol)
valid.genes <- c(valid.genes, unique(p2gGR$symbol)) %>% unique()

noLinks <- valid.genes[valid.genes %ni% names(p2gFreqs)]
zilch <- rep(0, length(noLinks))
names(zilch) <- noLinks
p2gFreqs <- c(p2gFreqs, zilch)
x <- 1:length(p2gFreqs)
rank_df <- data.frame(npeaks=p2gFreqs, rank=x)

# Plot barplot of how many linked peaks per gene
thresh <- 30
threshNpeaks <- rank_df$npeaks
threshNpeaks[threshNpeaks>thresh] <- thresh
nLinkedPeaks <- getFreqs(threshNpeaks)

df <- data.frame(nLinkedPeaks=as.integer(names(nLinkedPeaks)), nGenes=nLinkedPeaks)

pdf(paste0(plotDir, "/nGenes_with_nLinkedPeaks.pdf"), width=8, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=median(rank_df$npeaks), linetype="dashed")
dev.off()

# Get distribution of gene to peak linkages 
valid.peaks <- allPeaksGR$peakName
g2pFreqs <- getFreqs(p2gGR$peakName)

noLinks <- valid.peaks[valid.peaks %ni% names(g2pFreqs)]
zilch <- rep(0, length(noLinks))
names(zilch) <- noLinks
g2pFreqs <- c(g2pFreqs, zilch)
x <- 1:length(g2pFreqs)
g2p_rank_df <- data.frame(ngenes=g2pFreqs, rank=x)

# Plot barplot of how many linked peaks per gene
thresh <- 5
threshNgenes <- g2p_rank_df$ngenes
threshNgenes[threshNgenes>thresh] <- thresh
nLinkedGenes <- getFreqs(threshNgenes)

df <- data.frame(nLinkedGenes=as.integer(names(nLinkedGenes)), nPeaks=nLinkedGenes)

pdf(paste0(plotDir, "/nPeaks_with_nLinkedGenes.pdf"), width=5, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=median(g2p_rank_df$npeaks), linetype="dashed")
dev.off()

##########################################################################################
# Stacked bar chart of types of peaks in p2g links
##########################################################################################

all_peak_types <- getFreqs(allPeaksGR$peakType)
linked_peak_types <- getFreqs(allPeaksGR[allPeaksGR$peakName %in% p2gGR$peakName]$peakType)
multi_peaks <- rownames(g2p_rank_df[g2p_rank_df$ngenes > 1,])
multi_peak_types <- getFreqs(allPeaksGR[allPeaksGR$peakName %in% multi_peaks]$peakType)

plot_df <- data.frame(all_peak_types, linked_peak_types, multi_peak_types)
plot_df <- apply(plot_df, 2, function(x) x/sum(x))
melt_df <- reshape2::melt(plot_df)
colnames(melt_df) <- c("peak_type", "peak_group", "value")

pdf(paste0(plotDir, "/peak_types_barplot.pdf"), width=5, height=6)
stackedBarPlot(melt_df, xcol=2, fillcol=1, ycol=3, 
               cmap=getColorMap(cmaps_BOR$sambaNight, n=4, type="quantitative"), barwidth=0.9)
dev.off()

###################################################################################################
# Plot Tracks of ALL peak to gene links for select super-enhancer genes
###################################################################################################

source(paste0(scriptPath, "/cluster_labels.R"))

clustOrder <- fineOrder %>% unname() %>% unlist()

# Identify peaks present in superenhancer set
young_se_files <- list.files(
  path="/oak/stanford/groups/wjg/skim/projects/LDA/resources/superEnhancers/Hnisz2013Cell_SEs",
  pattern="*.csv$",
  full.names = TRUE
)
cell_source <- str_replace(basename(young_se_files), "\\.csv$", "")
names(young_se_files) <- cell_source

young_se_dt <- lapply(names(young_se_files), function(x){
  dt <- fread(young_se_files[x], header=FALSE, sep=",", skip=1)
  dt$source <- x
  dt
}) %>% rbindlist()
colnames(young_se_dt) <- c("enhID", "chr", "start", "end", "refseq", "enhRank", "SupEnh", "H3K27acDens", "rpmperbp", "source")
young_se_dt <- young_se_dt[SupEnh == 1] # Restrict to only super enhancers

# Convert refseq IDs to gene symbols
library(org.Hs.eg.db)
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq_to_symbol <- select(
  org.Hs.eg.db, 
  keys=unique(young_se_dt$refseq), 
  columns=c("REFSEQ", "SYMBOL"), 
  keytype="REFSEQ"
)
ref_to_sym <- refseq_to_symbol$SYMBOL
names(ref_to_sym) <- refseq_to_symbol$REFSEQ
young_se_dt$symbol <- ref_to_sym[young_se_dt$refseq]
young_se_dt <- young_se_dt[!is.na(young_se_dt$symbol)]

# Get the top N super enhancers from each source:
topN <- 250

se_sources <- unique(young_se_dt$source)
top_SEs_young <- lapply(se_sources, function(s){
  young_se_dt[source == s] %>% arrange(enhRank) %>% slice(1:topN) %>% pull(symbol)
})
names(top_SEs_young) <- se_sources

top_SEs <- top_SEs_young

all_top_SEs <- unlist(top_SEs) %>% unique()

# Label genes as being SE-associated or not
rank_df$SE <- ifelse(rownames(rank_df) %in% all_top_SEs, 1, 0)

# Plot the first 100 HRGs
rank_df_top <- rank_df[1:100,]
top_genes <- rownames(rank_df_top)

# Tracks of genes:
p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = top_genes, 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500,
  minCells=100
)

saveRDS(p, paste0(plotDir, "/Tracks-Top100_HRGs.rds"))

plotPDF(plotList = p, 
        name = "Tracks-Top100_HRGs.pdf", 
        ArchRProj = atac_proj, 
        addDOC = FALSE, 
        width = 8, height = 20)

# Plot the first 100 HRGs that are superenhancers
rank_df_SE <- filter(rank_df, SE == 1)
superenhancer_genes <- rownames(rank_df_SE[1:100,])

# Tracks of genes:
p1 <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = superenhancer_genes, 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500,
  minCells=100
)

saveRDS(p1, paste0(plotDir, "/Tracks-Top100_SE_HRGs.rds"))

plotPDF(plotList = p1, 
        name = "Tracks-Top100_SE_HRGs.pdf", 
        ArchRProj = atac_proj, 
        addDOC = FALSE, 
        width = 8, height = 20)

# High threshold marker genes for each cluster
Idents(rna_proj) <- "FineNamedClust"
markers <- FindAllMarkers(rna_proj, assay = "RNA", logfc.threshold = 1.0, min.pct = 0.3, min.diff.pct = 0.6, only.pos = T, max.cells.per.ident = 200)
markers$diff.pct <- markers$pct.1 - markers$pct.2

saveRDS(markers, paste0(plotDir, "/markergenes.rds"))

cts <- unique(rna_proj$FineNamedClust) 
cts <- cts[cts %in% markers$cluster] %>% sort()

markergenes <- list()
for (ct in cts) {
  genes <- filter(markers, cluster == ct) %>% rownames()
  genes <- genes[1:5]
  genes <- genes[!is.na(genes)]
  markergenes <- c(markergenes, genes)
}

markergenes <- unlist(markergenes) %>% unique()

# Tracks of genes:
p2 <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = markergenes, 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500,
  minCells=100, 
)

saveRDS(p2, paste0(plotDir, "/Tracks-Top5_markergenes.rds"))

plotPDF(plotList = p2, 
        name = "Tracks-Top5_markergenes.pdf", 
        ArchRProj = atac_proj, 
        addDOC = FALSE, 
        width = 8, height = 20)

###################################################################################################
# Assess if 'highly regulated genes' are enriched for super enhancer linked genes
###################################################################################################

# Hnisz super enhancers 2013
young_se_files <- list.files(
  path="/oak/stanford/groups/wjg/skim/projects/LDA/resources/superEnhancers/Hnisz2013Cell_SEs",
  pattern="*.csv$",
  full.names = TRUE
)
cell_source <- str_replace(basename(young_se_files), "\\.csv$", "")
names(young_se_files) <- cell_source

young_se_dt <- lapply(names(young_se_files), function(x){
  dt <- fread(young_se_files[x], header=FALSE, sep=",", skip=1)
  dt$source <- x
  dt
}) %>% rbindlist()
colnames(young_se_dt) <- c("enhID", "chr", "start", "end", "refseq", "enhRank", "SupEnh", "H3K27acDens", "rpmperbp", "source")
young_se_dt <- young_se_dt[SupEnh == 1] # Restrict to only super enhancers

# Convert refseq IDs to gene symbols
library(org.Hs.eg.db)
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq_to_symbol <- select(
  org.Hs.eg.db, 
  keys=unique(young_se_dt$refseq), 
  columns=c("REFSEQ", "SYMBOL"), 
  keytype="REFSEQ"
)
ref_to_sym <- refseq_to_symbol$SYMBOL
names(ref_to_sym) <- refseq_to_symbol$REFSEQ
young_se_dt$symbol <- ref_to_sym[young_se_dt$refseq]
young_se_dt <- young_se_dt[!is.na(young_se_dt$symbol)]

# Get the top N super enhancers from each source:
topN <- 250

se_sources <- unique(young_se_dt$source)
top_SEs_young <- lapply(se_sources, function(s){
  young_se_dt[source == s] %>% arrange(enhRank) %>% slice(1:topN) %>% pull(symbol)
})
names(top_SEs_young) <- se_sources

top_SEs <- top_SEs_young

all_top_SEs <- unlist(top_SEs) %>% unique()

# Label genes as being SE-associated or not
rank_df$SE <- ifelse(rownames(rank_df) %in% all_top_SEs, 1, 0)

# Cutoff for defining highly regulated genes (HRGs)
cutoff <- 20

# Hypergeometric enrichment of SE-associated genes in highly-regulated genes
q <- sum(rank_df[rank_df$npeaks > cutoff,]$SE)  # q = number of white balls drawn without replacement
m <- length(all_top_SEs)                        # m = number of white balls in urn
n <- nrow(rank_df) - m                          # n = number of black balls in urn
k <- nrow(rank_df[rank_df$npeaks > cutoff,])    # k = number of balls drawn from urn
SEmLog10pval <- -phyper(q, m, n, k, lower.tail=FALSE, log.p=TRUE)/log(10)

# Get p-value for all sources:
all_sources <- names(top_SEs)
pvals <- sapply(all_sources, function(s){
  rank_df$SE <- ifelse(rownames(rank_df) %in% top_SEs[[s]], 1, 0)
  # Hypergeometric enrichment of SE-associated genes in highly-regulated genes
  q <- sum(rank_df[rank_df$npeaks > cutoff,]$SE)  
  m <- length(top_SEs[[s]])                       
  n <- nrow(rank_df) - m                          
  k <- nrow(rank_df[rank_df$npeaks > cutoff,])    
  -phyper(q, m, n, k, lower.tail=FALSE, log.p=TRUE)/log(10)
}) %>% sort() %>% rev()

# Plot enrichment of SE's from different sources:
plot_df <- data.frame(rank=1:length(pvals), mlog10padj=-log10(p.adjust(10**-pvals, method="fdr")))
rownames(plot_df) <- names(pvals)
to_label <- c("CD14", "CD3", "UCSD_Lung", "BI_Brain_Hippocampus_Middle", "HeLa", "UCSD_Pancreas", "CD34_fetal", "HepG2")
plot_df$label <- ifelse(rownames(plot_df) %in% to_label, rownames(plot_df), "")
plot_df$color <- "royalblue1"

p <- (ggplot(data=plot_df, aes(x=rank, y=mlog10padj))
      + ggrepel::geom_text_repel(
        data = plot_df[plot_df$label != "",], aes(x=rank, y=mlog10padj, label=label), 
        size = 3,
        nudge_x = 2,
        hjust = "outward",
        segment.size = 0.1,
        box.padding=0.5,
        min.segment.length = 0, # draw all segments
        max.overlaps = Inf, # draw all labels
        color = "black")
      + geom_point(aes(color=color), size=2)
      + scale_color_manual(values=c("royalblue1"))
      + theme_BOR(border=FALSE)
      + ylab("Hypergeometric Enrichment -log10(Padj)")
      + xlab("")
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio = 1.0,
              legend.position = "none", # Remove legend
              axis.text.x = element_text(angle = 90, hjust = 1))
)

pdf(paste0(plotDir, "/SE_enrichment_by_source.pdf"))
p
dev.off()

# Rank plot:
rank_df$color <- ifelse(rank_df$npeaks >= cutoff, "royalblue1", "black")

# Label select super enhancer - associated genes
label_genes <- c(
  "ETS1", "EGFL7", "EPAS1", "COL4A1", "ETS1", "NOTCH1", "TBX2", "TBX4",
  "COL4A2", "FOXF1", "FOXA2", "TBX5", "IKZF1", "RUNX3", "NKX2-1", "FGFR2"
) %>% unique()

rank_df$label <- ifelse(rownames(rank_df) %in% label_genes, rownames(rank_df), "")

p <- (ggplot(data=rank_df, aes(x=rank, y=npeaks))
      + ggrepel::geom_text_repel(
        data = rank_df[rank_df$label != "",], aes(x=rank, y=npeaks, label=label), 
        size = 3,
        nudge_x = 2,
        #direction = "x",
        hjust = "outward",
        segment.size = 0.1,
        box.padding=0.5,
        min.segment.length = 0, # draw all segments
        max.overlaps = Inf, # draw all labels
        color = "black")
      + geom_point_rast(aes(color=color), size = 1)
      + scale_color_manual(values=c("black", "royalblue1"))
      + theme_BOR(border=FALSE)
      + ylab("N linked peaks")
      + xlab("")
      + ggtitle(sprintf("SE gene enrichment in top %s genes \n -log10(p-value) = %s", k, round(SEmLog10pval,2)))
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio = 1.0,
              legend.position = "none", # Remove legend
              axis.text.x = element_text(angle = 90, hjust = 1)) 
)

pdf(paste0(plotDir, "/nLinkedPeaksPerGene_rastr.pdf"))
p
dev.off()

saveRDS(rank_df, paste0(wd, "/HRG_rank_df.rds"))

#################################################################
# GO term enrichment for HRGs
#################################################################
source(paste0(scriptPath, "/GO_wrappers.R"))

rank_df <- readRDS(paste0(wd, "/HRG_rank_df.rds"))
cutoff <- 20
rank_df$hrg <- ifelse(rank_df$npeaks >= cutoff, TRUE, FALSE)

all_genes <- rownames(rank_df) %>% sort()
#hrgs <- rank_df[rank_df$hrg,] %>% rownames()
hrgs <- rank_df[1:50, ] %>% rownames()

upGO_MF <- calcTopGo(all_genes, interestingGenes=hrgs, nodeSize=5, ontology="MF") 
upGO_BP <- calcTopGo(all_genes, interestingGenes=hrgs, nodeSize=5, ontology="BP") 
upGO_CC <- calcTopGo(all_genes, interestingGenes=hrgs, nodeSize=5, ontology="CC") 

upGO_BP <- upGO_BP[order(as.numeric(upGO_BP$pvalue), decreasing = FALSE),]
upGO_MF <- upGO_MF[order(as.numeric(upGO_MF$pvalue), decreasing = FALSE),]
upGO_CC <- upGO_CC[order(as.numeric(upGO_CC$pvalue), decreasing = FALSE),]

pdf(paste0(plotDir, "/hrg_top100_GO_MFonly.pdf"), width = 8, height = 3)
topGObarPlot(upGO_MF, cmap = cmaps_BOR$comet, 
             nterms=3, border_color="black", 
             barwidth=0.85, barLimits=c(0, 6))
dev.off()

##########################################################################################
# Plot cell type specific P2G examples
##########################################################################################

clustOrder <- fineOrder$Stromal %>% unlist() %>% unname()
clustOrder <- clustOrder[which(clustOrder %ni% c("Schw", "Meso"))]

p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  features = flist,
  pal = FineNamedClustCmap,
  geneSymbol = c("ACAN"), 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500
)

pdf(paste0(plotDir, "/BrowserTrack_ACAN.pdf"), w = 6, h = 6)
grid::grid.newpage()
grid::grid.draw(p$ACAN)
dev.off()

#clustOrder <- fineOrder$Immune %>% unlist() %>% unname()
#clustOrder <- clustOrder[which(clustOrder %ni% c("Schw", "Meso"))]


# Violin plot for RNA expression
clustOrder <- fineOrder$Stromal %>% unlist() %>% unname()
clustOrder <- clustOrder[which(clustOrder %ni% c("Schw", "Meso"))]
clustsToPlot <- clustOrder
label_genes <- c("ACAN", "CCL3")
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

pList <- list()
for(gn in label_genes){
  gn <- "ACAN"
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
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Violin_Expression_ACAN.pdf"), width=6, height=2.5)
pList
dev.off()

# Overlap with known enhancers from Hu et al. (2012)
acan_p2g <- getP2G_GR(atac_proj, corrCutoff=0.45, filtNA = TRUE)
acan_p2g <- acan_p2g[which(acan_p2g$symbol == "ACAN")]

create_granges_list <- function(coord_strings) {
  # Initialize empty vectors to store the results
  chroms <- c()
  starts <- c()
  ends <- c()
  
  # Loop through each coordinate string
  for (coord_string in coord_strings) {
    matches <- regmatches(coord_string, regexec("^(chr[0-9XYM]+):(\\d+)-(\\d+)$", coord_string))
    
    if (length(matches[[1]]) == 4) {
      chroms <- c(chroms, matches[[1]][2])
      starts <- c(starts, as.numeric(matches[[1]][3]))
      ends <- c(ends, as.numeric(matches[[1]][4]))
    } else {
      stop(paste("Invalid coordinate format:", coord_string))
    }
  }
  
  # Create a GRanges object
  gr <- GRanges(seqnames = chroms, ranges = IRanges(start = starts, end = ends))
  return(gr)
}

hu_enhancers <- readr::read_csv("./hu_enhancers_hg38.csv", col_names = "coord")
hu_enhancers <- create_granges_list(hu_enhancers$coord)

# Overlap between p2gs and enhancers
ol <- findOverlaps(hu_enhancers, acan_p2g)
hu_enhancers_ol <- hu_enhancers[queryHits(ol)]
acan_p2g_ol <- acan_p2g[subjectHits(ol)]
  
enh_acan <- readr::read_csv("./hu_enhancers.csv")
mcols(hu_enhancers) <- enh_acan[3:6]
# Plot browser track with the enhancers and overlaps
flist <- list()
flist[["peaks"]] <- getPeakSet(atac_proj)
#flist[["Enhancers from Hu et al."]] <- hu_enhancers %>% resize(250, fix="center")
#flist[["Hu_enhancers"]] <- acan_p2g_ol
flist[["Hu_enhancers"]] <- hu_enhancers_ol
#flist[["linked_SNPs"]] <- unique(dis_full_fm_gr[dis_full_fm_gr$index_SNP %in% index_map]) %>% resize(250, fix="center")
#flist[["candidate_SNPs"]] <- candidate_GR[!duplicated(candidate_GR)] %>% resize(250, fix="center")

browser_region <- GRanges(
  seqnames = "chr15",
  ranges = IRanges(start = 88753442, end = 88853443)
)
ol1 <- findOverlaps(hu_enhancers, browser_region)
hu_enhancers_browser <- hu_enhancers[queryHits(ol1)]
hu_enhancers_confirmed <- hu_enhancers_browser[!is.na(hu_enhancers_browser$DevelopmentActiveEnhancer),]

hu_enhancers_confirmed <- hu_enhancers_confirmed[!is.na(hu_enhancers_confirmed$DevelopmentOnly),]
#flist[["Hu_confirmed"]] <- hu_enhancers_confirmed
flist[["Hu_confirmed"]] <- hu_enhancers_confirmed

p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  features = flist,
  pal = FineNamedClustCmap,
  geneSymbol = c("ACAN"), 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500
)

pdf(paste0(plotDir, "/BrowserTrack_ACAN_final.pdf"), w = 6, h = 6)
grid::grid.newpage()
grid::grid.draw(p$ACAN)
dev.off()

##########################################################################################
# Plot Peak2Gene heatmap
##########################################################################################

################################
nclust <- 30
p <- plotPeak2GeneHeatmap(
  atac_proj, 
  corCutOff = corrCutoff, 
  groupBy="FineNamedClust", 
  nPlot = 1000000, returnMatrices=FALSE, 
  k=nclust, seed=1, palGroup=FineNamedClustCmap
)
pdf(paste0(plotDir, sprintf("/peakToGeneHeatmap_LabelClust_k%s.pdf", nclust)), width=16, height=10)
p
dev.off()
################################

# Need to force it to plot all peaks if you want to match the labeling when you 'returnMatrices'.
p2gMat <- plotPeak2GeneHeatmap(
  atac_proj, 
  corCutOff = corrCutoff, 
  groupBy="FineNamedClust",
  nPlot = 1000000, returnMatrices=TRUE, 
  k=nclust, seed=1)

# Get association of peaks to clusters
kclust_df <- data.frame(kclust=p2gMat$ATAC$kmeansId, peakName=p2gMat$Peak2GeneLinks$peak, gene=p2gMat$Peak2GeneLinks$gene)

# Fix peakname
kclust_df$peakName <- sapply(kclust_df$peakName, function(x) strsplit(x, ":|-")[[1]] %>% paste(.,collapse="_"))

# Get motif matches
matches <- getMatches(atac_proj, "Motif")
r1 <- SummarizedExperiment::rowRanges(matches)
rownames(matches) <- paste(seqnames(r1),start(r1),end(r1),sep="_")
matches <- matches[names(allPeaksGR)]

clusters <- unique(kclust_df$kclust) %>% sort()

enrichList <- lapply(clusters, function(x){
  cPeaks <- kclust_df[kclust_df$kclust == x,]$peakName %>% unique()
  ArchR:::.computeEnrichment(matches, which(names(allPeaksGR) %in% cPeaks), seq_len(nrow(matches)))
}) %>% SimpleList
names(enrichList) <- clusters

# Format output to match ArchR's enrichment output
assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
  d <- lapply(seq_along(enrichList), function(y){
    enrichList[[y]][colnames(matches),x,drop=FALSE]
  }) %>% Reduce("cbind",.)
  colnames(d) <- names(enrichList)
  d
}) %>% SimpleList
names(assays) <- colnames(enrichList[[1]])
assays <- rev(assays)
res <- SummarizedExperiment::SummarizedExperiment(assays=assays)

formatEnrichMat <- function(mat, topN, minSig, clustCols=TRUE){
  plotFactors <- lapply(colnames(mat), function(x){
    ord <- mat[order(mat[,x], decreasing=TRUE),]
    ord <- ord[ord[,x]>minSig,]
    rownames(head(ord, n=topN))
  }) %>% do.call(c,.) %>% unique()
  pMat <- mat[plotFactors,]
  prettyOrderMat(pMat, clusterCols=clustCols)$mat
}

pMat <- formatEnrichMat(assays(res)$mlog10Padj, 5, 10, clustCols=FALSE)
# Save maximum enrichment
tfs <- strsplit(rownames(pMat), "_") %>% sapply(., `[`, 1)
rownames(pMat) <- paste0(tfs, " (", apply(pMat, 1, function(x) floor(max(x))), ")")

pMat <- apply(pMat, 1, function(x) x/max(x)) %>% t()

pdf(paste0(plotDir, sprintf("/enrichedMotifs_kclust_p2gHM_k%s.pdf", nclust)), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  pMat, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  #top_annotation = ta,
  row_names_side = "left",
  width = ncol(pMat)*unit(0.5, "cm"),
  height = nrow(pMat)*unit(0.33, "cm"),
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]"
)
draw(hm)
dev.off()

# GO enrichments of top N genes per cluster 
# ("Top" genes are defined as having the most peak-to-gene links)
source(paste0(scriptPath, "/GO_wrappers.R"))

kclust <- unique(kclust_df$kclust) %>% sort()
all_genes <- kclust_df$gene %>% unique() %>% sort()

# Save table of top linked genes per kclust
nGOgenes <- 200
topKclustGenes <- lapply(kclust, function(k){
  kclust_df[kclust_df$kclust == k,]$gene %>% getFreqs() %>% head(nGOgenes) %>% names()
}) %>% do.call(cbind,.)
outfile <- paste0(plotDir, sprintf("/topN_genes_kclust_k%s.tsv", nclust))
write.table(topKclustGenes, file=outfile, quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)

GOresults <- lapply(kclust, function(k){
  message(sprintf("Running GO enrichments on k cluster %s...", k))
  clust_genes <- topKclustGenes[,k]
  upGO <- rbind(
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="BP") 
    #calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=upGenes, nodeSize=5, ontology="CC")
  )
  upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
})

names(GOresults) <- paste0("cluster_", kclust)

# Plots of GO term enrichments:
pdf(paste0(plotDir, sprintf("/kclust_GO_3termsBPonlyBarLim_k%s.pdf", nclust)), width=10, height=2.5)
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
# Marker peaks in HRGs
##########################################################################################
saveArchRProject(atac_proj, outputDirectory = paste0(wd, "/lda_v2_HRP"), load = FALSE)

hrp_atac_proj <- loadArchRProject(paste0(wd, "/lda_v2_HRP"))

rank_df <- readRDS(paste0(wd, "/HRG_rank_df.rds"))
hrgs <- rownames(rank_df[which(rank_df$npeaks >= 20),])

require(EnsDb.Hsapiens.v86)
matchedGeneSymbols <- select(EnsDb.Hsapiens.v86,
                             keys = hrgs,
                             columns = c("SYMBOL", "GENEBIOTYPE"),
                             keytype = "SYMBOL")
genesToKeep <- matchedGeneSymbols %>% dplyr::filter(GENEBIOTYPE == "protein_coding") %>% dplyr::select(SYMBOL) %>% unlist() %>% unname()

hrgs <- hrgs[which(hrgs %in% genesToKeep)] #1052 genes
saveRDS(hrgs, paste0(wd, "/HRG_rank_df_biotype_filtered.rds"))

mkpeaks <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_v2/Plots/lda_markerPeaks_FineNamedClust.rds")
markerList <- getMarkers(mkpeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

# Sort and filter P2G and remove redundant P2G
full_p2gGR <- readRDS(file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]

filt_p2gGR <- full_p2gGR[full_p2gGR$Correlation > corrCutoff & 
                           full_p2gGR$VarQATAC > varCutoffATAC & 
                           full_p2gGR$VarQRNA > varCutoffRNA]
filt_p2gGR <- filt_p2gGR[!duplicated(paste0(filt_p2gGR$peakName, "_", filt_p2gGR$symbol))] %>% sort()
filt_p2gGR <- filt_p2gGR[!is.na(filt_p2gGR$idxATAC)]

saveRDS(filt_p2gGR, paste0(wd, "/multilevel_p2gGR_sorted_filtered.rds"))
filt_p2gGR <- readRDS(paste0(wd, "/multilevel_p2gGR_sorted_filtered.rds"))

# Organize peaks linked to each HRG
peak.list <- lapply(hrgs, function(hrg) {
  gr <- filt_p2gGR[which(filt_p2gGR$symbol == hrg),]
  mcols(gr) <- NULL
  names(ranges(gr)) <- NULL
  gr
})
peak.list <- GRangesList(peak.list)
names(peak.list) <- hrgs

#saveRDS(peak.list, paste0(wd, "/PeakList_HRG_biotype_filtered.rds"))

hrp_atac_proj <- addPeakAnnotations(
  ArchRProj = hrp_atac_proj, 
  regions = peak.list, 
  name = "HRP",
  force = TRUE
  )

enrichRegions <- peakAnnoEnrichment(
  seMarker = mkpeaks,
  ArchRProj = hrp_atac_proj,
  peakAnnotation = "HRP",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichRegions[,unname(unlist(fineOrder))], 
  n=5, 
  #clusterCols = FALSE, # This is currently bugged and does not work
  transpose=FALSE, 
  cutOff=log10pCut
)

#draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
#plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-FineNamedClust", width=8, height=12, ArchRProj=atac_proj, addDOC=FALSE)

# Heatmap of motif enrichments
plot_mat <- plotEnrichHeatmap(enrichRegions[,unname(unlist(fineOrder))], n=5, transpose=FALSE, 
                              cutOff=log10pCut, returnMatrix=TRUE)

plot_mat <- prettyOrderMat(plot_mat[,unname(unlist(fineOrder))], clusterCols=FALSE, cutOff=1)$mat

fontsize <- 8
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=unname(unlist(fineOrder)),col=list(atac_cluster=FineNamedClustCmap), 
                        show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))

hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = BuenColors::jdb_palette("brewer_red", type = "continuous"),
  top_annotation = ta,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = fontsize),
  column_names_gp = gpar(fontsize = fontsize),
  width = ncol(plot_mat)*unit(0.3, "cm"),
  height = nrow(plot_mat)*unit(0.3, "cm"),
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]",
  border_gp = gpar(col="black") # Add a black border to entire heatmap
)

pdf(paste0(plotDir, "/MarkerPeak-Enrichment-HRPs-FineNamedClust.pdf"), width=13, height=18)
draw(hm)
dev.off()

# Heatmap of HRP enrichment in cell types transposed ###
plot_mat <- plotEnrichHeatmap(enrichRegions[,unname(unlist(fineOrder))], n=5, transpose=TRUE, 
                              cutOff=log10pCut, returnMatrix=TRUE)

plot_mat <- prettyOrderMat(plot_mat[unname(unlist(fineOrder)),], clusterCols=FALSE, cutOff=1)$mat

fontsize <- 8
ht_opt$simple_anno_size <- unit(0.25, "cm")
#ta <- HeatmapAnnotation(HRGs=colnames(plot_mat),
#                        show_legend=c(HRGs=FALSE), show_annotation_name = c(HRGs=FALSE))

hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = BuenColors::jdb_palette("brewer_red", type = "continuous"),
  #top_annotation = ta,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = fontsize),
  column_names_side = "top",
  column_names_gp = gpar(fontsize = fontsize),
  width = ncol(plot_mat)*unit(0.3, "cm"),
  height = nrow(plot_mat)*unit(0.3, "cm"),
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]",
  border_gp = gpar(col="black") # Add a black border to entire heatmap
)

pdf(paste0(plotDir, "/MarkerPeak-Enrichment-HRPs-FineNamedClust_transposed.pdf"), width=20, height=10)
draw(hm)
dev.off()

# Add custom deviation matrix
hrp_atac_proj <- addBgdPeaks(hrp_atac_proj, force = TRUE)

addArchRThreads(threads = 1) 
hrp_atac_proj <- addDeviationsMatrix(
  ArchRProj = hrp_atac_proj, 
  peakAnnotation = "HRP",
  force = TRUE
)

plotVarDev <- getVarDeviations(hrp_atac_proj, plot = TRUE, name = "HRPMatrix")

pdf(paste0(plotDir, "/variable_HRP_scores.pdf"))
plotVarDev
dev.off()

markerHRP <- getFeatures(hrp_atac_proj, useMatrix = "HRPMatrix")
markerHRP <- sort(grep("z:", markerHRP, value = TRUE))

sorted_HRP <- paste0("z:",hrgs) # take the n peaks rank and rank HRPs based on n peaks
indices <- match(markerHRP, sorted_HRP)
markerHRP <- markerHRP[indices]

hrp_atac_proj <- addImputeWeights(hrp_atac_proj, reducedDims = "PCA")

hrpToPlot <- markerHRP[1:100]

p <- plotEmbedding(
  ArchRProj = hrp_atac_proj, 
  colorBy = "HRPMatrix", 
  name = "z:ACAN",
  embedding = "customUMAP"#,
  #imputeWeights = getImputeWeights(hrp_atac_proj)
)

plotPDF(p, name = "Variable-HRP-Scores-UMAPs-w-Imputation", width = 5, height = 5, ArchRProj = hrp_atac_proj, addDOC = FALSE)

saveArchRProject(hrp_atac_proj, outputDirectory = paste0(wd, "/lda_v2_HRP"), load = FALSE)

# Plot relationship between n_peaks associated with HRG and HRP to gene expression correlation #
corGEM_HM <- correlateMatrices(
  ArchRProj = hrp_atac_proj,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "HRPMatrix",
  reducedDims = "PCA"
)

corGEM_HM <- as.data.frame(corGEM_HM)

filt_rank_df <- rank_df %>% dplyr::filter(npeaks >= 20)
filt_rank_df$gene <- rownames(filt_rank_df)
filt_rank_df <- as.data.frame(filt_rank_df)

cor.df <- left_join(filt_rank_df, 
                    corGEM_HM, 
                    by = c("gene" = "GeneExpressionMatrix_matchName"))

# Clean up NAs
cor.df <- cor.df[which(!is.na(cor.df$GeneExpressionMatrix_name)),]

pdf(paste0(plotDir, "/npeaksBycorGEM_HM_scatter.pdf"))
ggplot(cor.df, aes(x = npeaks, y = cor)) +
  geom_point()
dev.off()

### Browser tracks

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

clustOrder <- fineOrder %>% unlist() %>% unname()

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  pal = FineNamedClustCmap,
  geneSymbol = c("CRISPLD2"), 
  loops = getPeak2GeneLinks(proj),
  tileSize=500
)
grid::grid.newpage()
grid::grid.draw(p$CRISPLD2)

# Tracks of genes: ####
# (Define plot region based on bracketing linked peaks)
promoterGR <- promoters(getGenes(atac_proj))

mPromoterGR <- promoterGR[promoterGR$symbol %in% label_genes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% label_genes]

# Restrict to only loops linking genes of interest
plotLoops <- getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff, resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# Bracket plot regions around SNPs
plotRegions <- lapply(label_genes, function(x){
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  gr <- GRanges(
    seqnames = seqnames(gr)[1],
    ranges = IRanges(start=lims[1], end=lims[2])
  )
  gr
}) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
                      width=width(plotRegions) + 0.05*width(plotRegions), 
                      fix="center")

label_genes <- c("")
# Tracks of genes:
p <- plotBrowserTrack(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  useGroups = clustOrder,
  pal = FineNamedClustCmap,
  plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
  sizes = c(7, 0.2, 1.25, 2.5),
  geneSymbol = "SOX2", 
  loops = getPeak2GeneLinks(atac_proj),
  tileSize=500,
  minCells=100, 
)

grid::grid.draw(p$SOX2)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Interesting-Genes.pdf", 
        ArchRProj = atac_proj, 
        addDOC = FALSE, 
        width = 8, height = 20)

p <- plotBrowserTrack(atac_proj, groupBy = "FineNamedClust", 
                      useGroups = clustOrder,
                      pal = FineNamedClustCmap,
                      geneSymbol = "PROX1",
                      loops = getPeak2GeneLinks(atac_proj)$Peak2GeneLinks)

grid::grid.draw(p$PROX1)
