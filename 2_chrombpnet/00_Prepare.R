#!/usr/bin/env Rscript

########################################
# Prepare files for training
########################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
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
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = TRUE, recursive = TRUE)
setwd(wd)

plotDir <- paste0(wd, "/plots")
dir.create(plotDir, showWarnings = TRUE, recursive = TRUE)

# Misc options
addArchRGenome("hg38")

# Set genome
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

# Minimum values for each cluster
min_reads <- 5e6 # number of reads
min_cells <- 140 # number of cells (originally 200)

# Number of peaks to use per cell type
peak_width <- 1001 # Peak width to use 
total_peaks <- 75000
nmarker_peaks <- 7500

# Block scientific notations
options(scipen = 999)

# Load ArchR Project
atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda")

##########################################################################################
# Prepare fragments files per cluster
##########################################################################################

# Export fragments file at cluster level defined by groupBy parameter
# *Note that this automatically creates a GroupFragments directory in the ArchR Project directory
# and creates a fragments tsv.gz file for each pseudobulked cluster
frags <- getGroupFragments(ArchRProj = atac_proj, groupBy = "FineNamedClust")

# Move fragment files to wd
frag.dir <- paste0(wd, "/GroupFragments")
file.rename(from = "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda/GroupFragments",
            to = paste0(wd, "/GroupFragments"))

frags <- list.files(frag.dir, full.names = T)

# Calculate fragments per cluster and check how many fall under cutoff
cellTypes <- names(table(atac_proj$FineNamedClust))

readsPerCluster <- data.frame(cluster = cellTypes, total_reads = NA)

# Count the number of fragments per cluster
for (i in seq_along(frags)) {
  data <- fread(cmd = paste0("zcat ", frags[i]))
  readsPerCluster$cluster[i] <- cellTypes[i]
  readsPerCluster$total_reads[i] <- nrow(data)
}

# reorder and save tsv
readsPerCluster <- readsPerCluster %>% arrange(desc(total_reads))
write_tsv(readsPerCluster, paste0(plotDir, "/fragmentsPerCluster.tsv"))

# Plot fragments per cluster
readsPerCluster <- read_tsv(paste0(plotDir, "/fragmentsPerCluster.tsv")) %>% as.data.frame()
nPassing <- length(which(readsPerCluster$total_reads > min_reads))

pdf(paste0(plotDir, "/fragmentsPerCluster.pdf"), w = 10, h = 4)
ggplot(readsPerCluster, 
       aes(
         x = reorder(cluster, -total_reads), 
         y = total_reads)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = min_reads) + 
  theme_ArchR() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(1, 10e10)
  ) +
  xlab("") +
  ylab("Fragments per cluster") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14)
    ) +
  ggtitle(paste0("Number of clusters with >", min_reads, " reads per cluster = ", nPassing, "/", length(cellTypes)))
dev.off()

##########################################################################################
# Prepare cluster specific peaks
##########################################################################################

# Get all peaks from ArchR project
peaks_gr <- getPeakSet(atac_proj)
peaks_gr$peakName <- (peaks_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
seqinfo(peaks_gr) <- seqinfo(genome)[seqlevels(peaks_gr)]

# Get blacklist regions used by chrombpnet
blacklist <- import.bed("/oak/stanford/groups/wjg/skim/projects/LDA/resources/chrombpnet/blacklist.bed.gz")

# Extend blacklist by 1057bp on both sides
blacklist_ext <- resize(blacklist, width = width(blacklist) + (1057 * 2), fix = "center")

# Get peaks that were originally called for each cluster
cellTypes <- names(table(atac_proj$FineNamedClust))

# Create directory to store peak files and the blacklist regions used to filter
peaksDir <- paste0(wd, "/peaks")
dir.create(peaksDir, showWarnings = FALSE, recursive = TRUE)
export(blacklist_ext, paste0(wd, "/blacklist_extended.bed"), "bed")

# For each cluster create a bed file of peaks called for that cluster
for (cluster in cellTypes) {
  # Extract peaks called for that cluster
  peaks <- getClusterPeaks(atac_proj, clusterNames=cluster, peakGR=peaks_gr, groupBy = "FineNamedClust")
  
  # Filter peaks that overlap with the extended blacklist sites
  peaks <- peaks[!peaks %over% blacklist_ext]
  
  # Create a bed file for that cluster specific peaks
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
    # and chrombpnet takese in 1000bp peak regions
    )
  message(paste0("writing bed file for ", cluster))
  write.table(bed, file = paste0(peaksDir, "/", cluster, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
}

# Count number of peaks per cluster identified
nPeaksPerCluster <- data.frame(cluster = cellTypes, nPeaks = NA)
for (i in seq_along(cellTypes)) {
  cluster <- cellTypes[i]
  bedfile <- paste0(peaksDir, "/", cluster, ".bed")
  nPeaks <- readLines(bedfile) %>% length()
  nPeaksPerCluster$cluster[i] <- cluster
  nPeaksPerCluster$nPeaks[i] <- nPeaks
}
nPeaksPerCluster <- nPeaksPerCluster %>% arrange(desc(nPeaks))
write_tsv(nPeaksPerCluster, paste0(plotDir, "/nPeaksPerCluster.tsv"))

# Plot nFrags by nPeaks for each cluster
perClusterMetadata <- full_join(nPeaksPerCluster, readsPerCluster, by = "cluster")
write_tsv(perClusterMetadata, paste0(plotDir, "/perClusterMetadata.tsv"))

pdf(paste0(plotDir, "/nFragsBynPeaksPerCluster.pdf"), w = 10, h = 6)
ggplot(perClusterMetadata, aes(x = total_reads, y = nPeaks)) + 
  geom_vline(xintercept = 5e6, linetype = "dashed") +
  annotate("text", x=5.2e6, y=40000, label="5e6 fragments", angle=90) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = cluster), max.overlaps = 15) +
  xlab("nFragments per Cluster") +
  ylab("nPeaks per Cluster") +
  scale_x_log10() + theme_bw()
dev.off()

##########################################################################################
# Prepare bigwigs for each cluster
##########################################################################################

atac_proj <- getGroupBW(
  ArchRProj = atac_proj,
  groupBy = "FineNamedClust",
  normMethod = "ReadsInTSS",
  tileSize = 500,
  maxCells = 1000,
  ceiling = 4
)












########### Testing #########

##########################################################################################
# Prepare marker peaks for each cluster (marker peaks + peaks called for that cluster)
##########################################################################################

# Get all peaks from ArchR project
peaks_gr <- getPeakSet(atac_proj)
peaks_gr$peakName <- (peaks_gr %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
seqinfo(peaks_gr) <- seqinfo(genome)[seqlevels(peaks_gr)]

# Identify blacklist regions in genome
use_chr <- unique(seqnames(peaks_gr))
use_masks <- c("AGAPS", "AMB")

# Get blacklist regions used by chrombpnet
blacklist <- import.bed("/oak/stanford/groups/wjg/skim/projects/LDA/resources/chrombpnet/blacklist.bed.gz")
# Extend blacklist by 1057bp on both sides
blacklist_ext <- resize(blacklist, width = width(blacklist) + (1057 * 2), fix = "center")
export(blacklist_ext, paste0(wd, "/blacklist_extended.bed"), "bed")

# Obtain full mask GR
masked_gr <- lapply(use_chr, function(chr){
  message(sprintf("Getting masks from %s...", chr))
  masks <- Biostrings::masks(genome[[as.character(chr)]])
  valid_masks <- masks@NAMES[sapply(masks@nir_list, length)>0]
  valid_masks <- valid_masks[valid_masks %in% use_masks]
  full_range <- sapply(valid_masks, function(msk){GRanges(chr, as(masks[[msk]], "IRanges"))}) %>%
    as(., "GRangesList") %>% unlist()
  full_range$mask <- names(full_range)
  full_range
  }) %>% as(., "GRangesList") %>% unlist() %>% sort()

# Create a union blacklist set
blacklist_combined <- reduce(union(blacklist, masked_gr))
export(blacklist_combined, paste0(wd, "/blacklist_combined.bed"), "bed")

# Filter cell types for number of cells per cluster
cell_types <- getFreqs(atac_proj$CustomNamedClust)
cell_types <- names(cell_types)[cell_types > min_cells]

# Identify marker peaks for each retained cluster
# We will use the 'top' N peaks (by score) for each cluster as well as all
# 'marker peaks' for that cluster to enhance cell-type specificity of model training data

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = atac_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "FineNamedClust",
    useGroups = cell_types,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff="FDR <= 0.1 & Log2FC >= 0.5")

plot_df <- lapply(names(markerList), function(name) {
  data.frame(Cluster = name, nMarkerPeaks = nrow(markerList[[name]]))
})

# visual check on distribution of number of marker peaks per cluster
plot_df <- bind_rows(plot_df)
plot_df <- plot_df %>% arrange(desc(nMarkerPeaks))
hist(plot_df$nMarkerPeaks)

# Select the top nmarker_peaks for each cluster
markerList <- lapply(markerList, function(mdf){
  mdf$peakName <- paste0(mdf$seqnames, "_", mdf$start, "_", mdf$end)
  head(mdf, nmarker_peaks) # Subset to most significant N marker peaks
  })

# Identify 'top' peaks for cell type and combine with significant set of marker peaks
cts_peaks <- lapply(cell_types, function(ct){
  message(sprintf("Gettting peaks from cluster %s...", ct))
  peaks <- getClusterPeaks(atac_proj, 
                           clusterNames=c(ct), 
                           peakGR=peaks_gr, 
                           originalScore=TRUE,
                           groupBy = "FineNamedClust")
  peaks <- peaks[!peaks %over% blacklist_combined]
  # Sort peaks by score
  peaks <- peaks[order(peaks$score, decreasing=TRUE)]
  # Prioritize 'marker' peaks for each cell type
  marker_peak_names <- markerList[[ct]]$peakName
  peaks <- c(peaks[peaks$peakName %in% marker_peak_names], peaks)
  # Deduplicate and clean up
  peaks <- peaks[!duplicated(peaks$peakName)]
  peaks <- resize(peaks,width=peak_width, fix="center") %>% trim_oob() %>% trim_N_seqs(.,genome)
  peaks <- peaks %>% head(total_peaks) %>% sort()
  peaks$GC <- gcContent(peaks, genome)
  peaks
  })

names(cts_peaks) <- cell_types

# Only keep cell types that have at least 48,750 peaks total
cts_peaks <- cts_peaks[unlist(lapply(cts_peaks, length)) > total_peaks*0.65]

# Plot similarity between input training data (Jaccard index)

# Grab cell types
cell_types <- names(cts_peaks)

# Function to calculate jaccard of peaks used in each model
jaccard <- function(v1, v2){
  # Calculate jaccard index between two vectors v1 and v2
  length(intersect(v1,v2))/length(unique(c(v1,v2)))
}

jacc_mat <- matrix(0, length(cell_types), length(cell_types))
for(i in 1:nrow(jacc_mat)){
  for(j in 1:ncol(jacc_mat)){
    v1 <- cts_peaks[[i]]$peakName
    v2 <- cts_peaks[[j]]$peakName
    jacc_mat[i,j] <- jaccard(v1, v2)
  }
}
rownames(jacc_mat) <- cell_types
colnames(jacc_mat) <- cell_types

jacc_mat[jacc_mat == 1.0] <- NA

pmat <- jacc_mat

# Add cluster labels and order the matrix
source(paste0(scriptPath, "/cluster_labels.R"))
clustOrder <- unlist(FineNamedClust) %>% unique()
clustOrder <- clustOrder[clustOrder %in% colnames(pmat)]
pmat <- pmat[clustOrder, clustOrder]

# Jaccard heatmap
pdf(paste0(plotDir, "/TrainingPeaks_jaccard_index_1000bp.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  pmat,
  limits=c(0,0.75),
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  column_title="FineNamedClust", 
  column_title_side="bottom",
  row_title="FineNamedClust",
  dataColors = cmaps_BOR$solar_extra,
  row_names_side = "left",
  width = ncol(pmat)*unit(0.5, "cm"),
  height = nrow(pmat)*unit(0.5, "cm"),
  legendTitle="Jaccard Index",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()

# Create directory to store peak files and the blacklist regions used to filter
peaksDir <- paste0(wd, "/peaks_markers")
dir.create(peaksDir, showWarnings = FALSE, recursive = TRUE)

# For each cluster create a bed file
for (cluster in names(cts_peaks)) {
  # Create a bed file for that cluster specific peaks
  peaks <- cts_peaks[[cluster]]
  bed <- data.frame(
    seqnames = seqnames(peaks),
    starts = start(peaks), #BED files are 0-indexed but granges are not
    ends = end(peaks),
    names = peaks$peakName,
    scores = peaks$score,
    strands = ".",
    col_6 = ".",
    col_7 = ".",
    col_8 = ".",
    summit = 500 #chrombpnet uses the 9th column as the summit position to calculate gc content. 500 since ArchR creates fixed width peaks that are 501bp long
    # and chrombpnet uses 1000bp peak regions
  )
  write.table(bed, file = paste0(peaksDir, "/", cluster, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
}








# Temp
p <- plotGroups(
  ArchRProj = atac_proj, 
  groupBy = "FineNamedClust", 
  colorBy = "cellColData", 
  name = "FRIP",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
