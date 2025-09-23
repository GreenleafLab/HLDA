#!/usr/bin/env Rscript

########################################
# Summarize Modisco outputs
########################################

# Load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(rvest) # for html file manipulation
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
#atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda")
#rna_obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

################################################
# Organize modisco outputs for each cluster
################################################
modisco_dir <- paste0(wd, "/modisco_outputs")

clusters <- sub("_reports$", "", basename(list.dirs(modisco_dir, recursive = F))) #Grab cluster names from directory names in modisco_dir

df.list <- list()

# For each modisco report, extract predicted motifs and qvalues
for (cluster in clusters) {
  report_html <- read_html(paste0(modisco_dir, "/", cluster, "_reports/motifs.html"))
  modisco_df <- report_html %>% html_node("table") %>% html_table()
  df <- dplyr::select(modisco_df, c("num_seqlets", "match0", "qval0", "match1", "qval1", "match2", "qval2"))
  df$cluster <- cluster
  df.list <- c(df.list, list(df))
}

# Combine into one df
modisco_df <- do.call(rbind, df.list)
saveRDS(modisco_df, paste0(plotDir, "/modisco_outputs_df.rds"))
#modisco_df <- readRDS(paste0(plotDir, "/modisco_outputs_df.rds"))
################################################
# Plot motif enrichment heatmap
################################################
#enrichMotifs <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/Figure01")

# Name adjustments for motif duplicates and mouse motif names
modisco_df$motif_name <- sub("^(.*?)_.*$", "\\1", modisco_df$match0)
modisco_df$motif_name <- str_to_upper(modisco_df$motif_name)
modisco_df$motif_name <- gsub("\\.MOUSE", "", modisco_df$motif_name)

modisco_df$log10p <- -log10(modisco_df$qval0) # -log10 of qval

# Subset to motifs that have at least some enrichment
log10p_cutoff <- 5
filt.modisco_df <- modisco_df %>% filter(log10p > log10p_cutoff) # note the comparison is flipped since -log10 of adj P value

# Transform the df to motif by cluster (in cases of duplicate motif names, choose the qval that's more significant)
filt.modisco_df <- filt.modisco_df %>% 
  dplyr::select("motif_name","log10p", "cluster") %>% 
  pivot_wider(id_cols = motif_name, names_from = cluster, values_from = log10p, values_fn = max, values_fill = 0) %>% as.data.frame()

# Clean up rowname
rownames(filt.modisco_df) <- filt.modisco_df$motif_name
filt.modisco_df <- filt.modisco_df[,-1]

# Reorder matrix
plot_mat <- as.matrix(filt.modisco_df)
fineOrder <- unname(fineOrder) %>% unlist()
clustOrder <- fineOrder[which(fineOrder %in% colnames(plot_mat))]
plot_mat <- prettyOrderMat(plot_mat[,clustOrder], clusterCols=FALSE, cutOff=1)$mat
plot_mat_df <- plot_mat %>% as.data.frame()

# Plot heatmap
pdf(paste0(plotDir, "/MarkerPeak-ChromBPNET-MotifEnriched-Heatmap-FineNamedClust.pdf"), width=13, height=18)
fontsize <- 8
ht_opt$simple_anno_size <- unit(0.25, "cm")
#ta <- HeatmapAnnotation(atac_cluster=clustOrder,col=list(atac_cluster=FineNamedClustCmap), 
                        show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat_df, 
  limits=c(0.0,20.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  #top_annotation = ta,
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

saveRDS(plot_mat_df, paste0(plotDir,"/MarkerPeak-ChromBPNET-MotifEnriched-Heatmap-FineNamedClust.rds"))


