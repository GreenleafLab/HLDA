#!/usr/bin/env Rscript

##########################################################################
# Link genes to gkmSVM model prioritized SNPs
##########################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(parallel)
})

# Set Threads to be used
ncores <- 12

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda"
gkm_res_dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM/snp_results"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
addArchRGenome("hg38")

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
plotDir <- paste0(gkm_res_dir, "/plots") #paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

# Color Maps
compartmentCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/lungClusterColors.rds") %>% unlist()
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds") %>% unlist()
#BroadNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_BroadNamedClust_cmap.rds") %>% unlist()

sample_cmap <- readRDS(paste0("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/", "sample_cmap.rds"))
#rna_sample_cmap <- sample_cmap[names(sample_cmap) %in% unique(rna_proj$Sample)] %>% unlist()

gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")
#gest_age_cmap <- gest_age_cmap[names(gest_age_cmap) %in% unique(rna_proj$age)]

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Correlation cutoff for identifying p2g linkages
corrCutoff <- 0.45

# Retrieve GEmat prior to re-assigning p2g links
GEmat <- getMatrixFromProject(atac_proj, useMatrix="GeneExpressionMatrix")

##########################################################################################
# Links fine-mapped SNPs to candidate genes using Peak-to-Gene links
##########################################################################################

fmGR <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM/snp_fastas/250bpSNPCentered.rds")

# Load full project p2g links, plot loops, etc.
full_p2gGR <- readRDS(file="/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/P2G_Analysis/multilevel_p2gGR.rds") # NOT merged or correlation filtered
plot_loop_list <- readRDS(file="/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/P2G_Analysis/multilevel_plot_loops.rds")

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$symbol, "-", full_p2gGR$peakName))] %>% sort()
filt_p2gGR <- filt_p2gGR[!is.na(filt_p2gGR$idxATAC)] # Remove ATAC peak indices that are NA for some reason

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

# Get full merged p2g links
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)

pol <- findOverlaps(p2gGR, fmGR, type="any", maxgap=-1L, ignore.strand=TRUE)
expandFMGR <- fmGR[to(pol)]
expandFMGR$linkedGene <- p2gGR[from(pol)]$symbol
expandFMGR$SNP_to_gene <- paste(expandFMGR$linked_SNP, expandFMGR$linkedGene, sep="_")

##########################################################################################
# Read in results of gkmSVM models
##########################################################################################

# Load intermediate results
full_sig_res <- readRDS(paste0(gkm_res_dir, "/significant_hits_table.rds"))
snp_gkmexplain <- readRDS(paste0(gkm_res_dir, "/snp_gkmexplain.rds"))

# Filter fine-mapping gene links by those that were in gkmSVM significant hits
keep_dis <- c("Asthma (childhood onset)", "Asthma", "Asthma (adult onset)", "Childhood asthma with severe exacerbations", "Lung disease severity in cystic fibrosis",
  "Lung function (FVC)", "FEV1", "Lung function (FEV1/FVC)", "FEV1FVC")
dis_FM_GR <- expandFMGR[expandFMGR$disease_trait %in% keep_dis]

# Add linked genes to gkmSVM results
snp_to_gene_df <- mcols(expandFMGR) %>% as.data.frame() %>% group_by(linked_SNP) %>% summarize(genes=paste(unique(linkedGene), collapse=";")) %>% as.data.frame()
gene_vec <- snp_to_gene_df$genes
names(gene_vec) <- snp_to_gene_df$linked_SNP

full_sig_res$linkedGenes <- gene_vec[full_sig_res$snp]

readr::write_tsv(full_sig_res, file = "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM/snp_results/full_sig_res.tsv")

##########################################################################################
# Plot seqlogos using original gkmexplain matrix
##########################################################################################

library(ggseqlogo)

plotSNPimportance <- function(snp, snp_table, gkm_explain_output, celltype, gr, indices=101:151){
  # Plot the ref, alt, and delta importance scores for a given SNP in a given cell type
  # snp = the snp to plot
  # snp_table = df of snp hits with the following columns: (region, snp, ref_score, alt_score, score_delta,...)
  # gkm_explain_output = giant list of full gkmexplain output 
  # celltype = celltype to plot
  # gr = genomic range of snps

  # Get required names for accessing data
  snp_info <- snp_table[(snp_table$snp == snp & snp_table$cluster == celltype),]
  ref_group <- paste0(celltype, "-ref_snp_seqs")
  alt_group <- paste0(celltype, "-alt_snp_seqs")
  region <- snp_info$region[1]
  region <- paste0(region, "_", snp, "_125") # Need to adjust this for different relative SNP position

  # Get SNP region
  snp_gr <- gr[gr$linked_SNP == snp] %>% resize(width=50, fix="center")
  true_region <- (snp_gr %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})
  ref_alt <- snp_gr$linked_refalt[1]
  disease <- snp_gr$disease_trait[1]

  # Get info on suspected motif match, linked genes, etc
  linked_genes <- snp_info$linkedGenes[1]
  top_motifs <- snp_info$top_motifs[1]
  fm_prob <- snp_info$fm_probs[1]

  # Get matrices for plotting
  ref_matrix <- gkm_explain_output$gkmexplain_output[[ref_group]]$seq_matrices[[region]][,indices]
  alt_matrix <- gkm_explain_output$gkmexplain_output[[alt_group]]$seq_matrices[[region]][,indices]
  delta_matrix <- ref_matrix - alt_matrix

  # Generate plot for all matrices
  upper_lim <- max(rbind(ref_matrix, alt_matrix, delta_matrix)) * 1.2
  lower_lim <- min(rbind(ref_matrix, alt_matrix, delta_matrix)) * 1.2
  mat_list <- list("ref_allele"=ref_matrix, "alt_allele"=alt_matrix, "delta"=delta_matrix)

  p <- (
    ggseqlogo(mat_list, method='custom', seq_type='dna', ncol=1) 
    + ylab('gkmexplain importance')
    + ylim(lower_lim, upper_lim)
    + ggtitle(paste0(disease, " - ", celltype, " - ", snp, " | ", ref_alt, " - ", true_region, "\n", 
      linked_genes, "\n", top_motifs, "\n",  "FM probabiltiy: ", fm_prob))
    )
  p
}

valid_asthma_snps <- fmGR[fmGR$disease_trait %in% c("Asthma (childhood onset)", "Asthma", "Asthma (adult onset)", "Childhood asthma with severe exacerbations")]$linked_SNP
valid_spirometry_snps <- fmGR[fmGR$disease_trait %in% c("Lung function (FVC)", "FEV1", "Lung function (FEV1/FVC)", "FEV1FVC")]$linked_SNP
valid_cf_snps <- fmGR[fmGR$disease_trait %in% c("Lung disease severity in cystic fibrosis")]$linked_SNP

asthma_sig_res <- full_sig_res[full_sig_res$snp %in% valid_asthma_snps,]
asthma_sig_res <- asthma_sig_res[order(abs(asthma_sig_res$prominence), decreasing=TRUE),]

spirometry_sig_res <- full_sig_res[full_sig_res$snp %in% valid_spirometry_snps,]
spirometry_sig_res <- spirometry_sig_res[order(abs(spirometry_sig_res$prominence), decreasing=TRUE),]

cf_sig_res <- full_sig_res[full_sig_res$snp %in% valid_cf_snps,]
cf_sig_res <- cf_sig_res[order(abs(cf_sig_res$prominence), decreasing=TRUE),]

plot_asthma <- asthma_sig_res[!duplicated(asthma_sig_res$snp),]
plot_spirometry <- spirometry_sig_res[!duplicated(spirometry_sig_res$snp),]
plot_cf <- cf_sig_res[!duplicated(cf_sig_res$snp),]

# Only plot hits that have a linked gene
plot_asthma <- plot_asthma[!is.na(plot_asthma$linkedGenes),]
plot_spirometry <- plot_spirometry[!is.na(plot_spirometry$linkedGenes),]
plot_cf <- plot_cf[!is.na(plot_cf$linkedGenes),]

# Most prominent, high-effect SNP hits
nplot <- min(nrow(plot_asthma), 100)
pdf(paste0(gkm_res_dir, "/most_prominent_asthma_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nplot){
  snp <- plot_asthma$snp[i]
  ct <- plot_asthma$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, plot_asthma, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()

nplot <- min(nrow(plot_spirometry), 100)
pdf(paste0(gkm_res_dir, "/most_prominent_spirometry_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nplot){
  snp <- plot_spirometry$snp[i]
  ct <- plot_spirometry$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, plot_spirometry, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()

nplot <- min(nrow(plot_cf), 100)
pdf(paste0(gkm_res_dir, "/most_prominent_cf_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nplot){
  snp <- plot_cf$snp[i]
  ct <- plot_cf$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, plot_cf, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()

# Plotting all the potential candidate snp as tracks

all_candidate_snps <- c(plot_asthma$snp[1:100], plot_spirometry$snp[1:100]) %>% unique()

# Plot all model results for candidate SNPs:
sub_sig_res <- full_sig_res[sapply(full_sig_res$snp, function(x) grepl(paste(all_candidate_snps, collapse="|"), x)),]

# pdf(paste0(gkm_res_dir, "/all_candidate_snp_model_logos.pdf"), width=12, height=7)
# plotList <- list()
# for(i in 1:nrow(sub_sig_res)){
#   snp <- sub_sig_res$snp[i]
#   ct <- sub_sig_res$cluster[i]
#   plotList[[i]] <- plotSNPimportance(snp, sub_sig_res, snp_gkmexplain, celltype=ct, gr=fmGR)
# }
# plotList
# dev.off()


# Tracks of genes: 
disPicsGR <- expandFMGR[expandFMGR$disease_trait %in% keep_dis]
promoterGR <- promoters(getGenes(atac_proj))
candidate_GR <- disPicsGR[disPicsGR$linked_SNP %in% all_candidate_snps]
index_df <- data.frame(index=candidate_GR$index_SNP, candidate=candidate_GR$linked_SNP)
index_df <- index_df[!duplicated(index_df$candidate),]
index_map <- index_df$index
names(index_map) <- index_df$candidate

# The original index SNP is often not in a peak, so in order to find it we need to load the unfiltered fm gr
fm_dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/resources/gwas/PICS2"
full_fm_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))
dis_full_fm_gr <- full_fm_gr[full_fm_gr$disease_trait  %in% keep_dis]
index_GR <- dis_full_fm_gr[dis_full_fm_gr$linked_SNP %in% index_map] %>% unique()

# Marker Genes
markerGenes <- candidate_GR$linkedGene %>% unique()

# Combine plot loops from multi-level p2g linking
plotLoops <- unlist(as(plot_loop_list, "GRangesList"))
plotLoops <- plotLoops[order(plotLoops$value, decreasing=TRUE)] %>% unique() %>% sort()

mPromoterGR <- promoterGR[promoterGR$symbol %in% markerGenes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% markerGenes]

flist <- list()
flist[["peaks"]] <- getPeakSet(atac_proj)
flist[["index_SNPs"]] <- unique(index_GR) %>% resize(250, fix="center")
flist[["linked_SNPs"]] <- unique(dis_full_fm_gr[dis_full_fm_gr$index_SNP %in% index_map]) %>% resize(250, fix="center")
flist[["candidate_SNPs"]] <- candidate_GR[!duplicated(candidate_GR)] %>% resize(250, fix="center")

sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops <- plotLoops[width(plotLoops) > 500]

# Bracket plot regions around SNPs
plotRegions <- lapply(markerGenes, function(x){
  candidate_snps <- candidate_GR[candidate_GR$linkedGene == x]$linked_SNP
  gr <- c(
    range(candidate_GR[candidate_GR$linkedGene == x]), 
    range(index_GR[index_GR$linked_SNP %in% index_map[candidate_snps]]), 
    resize(range(mPromoterGR[mPromoterGR$symbol == x]), width=2000))
  lims <- grLims(gr)
  message(sprintf("Trying %s...", x))
  gr <- GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start=lims[1], end=lims[2])
    )
  gr
  }) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
  width=width(plotRegions) + 0.1*width(plotRegions), 
  fix="center")

plotRegions <- resize(plotRegions, width=ifelse(width(plotRegions) > 100000, width(plotRegions), 100000), fix="center")

FclustOrder <- fineOrder %>% unname() %>% unlist()

# FineNamedClust
p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "FineNamedClust",
    useGroups = FclustOrder,
    features = flist,
    pal = FineNamedClustCmap,
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(8, 0.8, 1.25, 0.5),
    region = plotRegions, 
    loops = getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff), # All peak-to-gene loops
    tileSize=250,
    minCells=100,
    title = markerGenes
)

plotPDF(plotList = p, 
    name = "candidate_snp_tracks_100k_width_allClust_allLoops.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 8, height = 10)


############# Specific snps only############
## Candidate SNPs of interest
all_candidate_snps <- c(
  "rs10975479", # IL33 linked to asthma
  "rs16970707", # HP linked to FVC
  "rs11067278", # TBX3 linked to FEV1/FVC
  "rs7578157" # ID2 linked to FEV1FVC
)

# Plot all model results for candidate SNPs:
sub_sig_res <- full_sig_res[sapply(full_sig_res$snp, function(x) grepl(paste(all_candidate_snps, collapse="|"), x)),]

pdf(paste0(gkm_res_dir, "/selected_candidate_snp_model_logos.pdf"), width=12, height=7)
plotList <- list()
for(i in 1:nrow(sub_sig_res)){
  snp <- sub_sig_res$snp[i]
  ct <- sub_sig_res$cluster[i]
  plotList[[i]] <- plotSNPimportance(snp, sub_sig_res, snp_gkmexplain, celltype=ct, gr=fmGR)
}
plotList
dev.off()


# Tracks of genes: 
disPicsGR <- expandFMGR[expandFMGR$disease_trait %in% keep_dis]
promoterGR <- promoters(getGenes(atac_proj))
candidate_GR <- disPicsGR[disPicsGR$linked_SNP %in% all_candidate_snps]
index_df <- data.frame(index=candidate_GR$index_SNP, candidate=candidate_GR$linked_SNP)
index_df <- index_df[!duplicated(index_df$candidate),]
index_map <- index_df$index
names(index_map) <- index_df$candidate

# The original index SNP is often not in a peak, so in order to find it we need to load the unfiltered fm gr
fm_dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/resources/gwas/PICS2"
full_fm_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))
dis_full_fm_gr <- full_fm_gr[full_fm_gr$disease_trait  %in% keep_dis]
index_GR <- dis_full_fm_gr[dis_full_fm_gr$linked_SNP %in% index_map] %>% unique()

# Marker Genes
markerGenes <- candidate_GR$linkedGene %>% unique()

# Combine plot loops from multi-level p2g linking
plotLoops <- unlist(as(plot_loop_list, "GRangesList"))
plotLoops <- plotLoops[order(plotLoops$value, decreasing=TRUE)] %>% unique() %>% sort()

mPromoterGR <- promoterGR[promoterGR$symbol %in% markerGenes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% markerGenes]

flist <- list()
flist[["peaks"]] <- getPeakSet(atac_proj)
flist[["index_SNPs"]] <- unique(index_GR) %>% resize(250, fix="center")
flist[["linked_SNPs"]] <- unique(dis_full_fm_gr[dis_full_fm_gr$index_SNP %in% index_map]) %>% resize(250, fix="center")
flist[["candidate_SNPs"]] <- candidate_GR[!duplicated(candidate_GR)] %>% resize(250, fix="center")

sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops <- plotLoops[width(plotLoops) > 500]

# Bracket plot regions around SNPs
plotRegions <- lapply(markerGenes, function(x){
  candidate_snps <- candidate_GR[candidate_GR$linkedGene == x]$linked_SNP
  gr <- c(
    range(candidate_GR[candidate_GR$linkedGene == x]), 
    range(index_GR[index_GR$linked_SNP %in% index_map[candidate_snps]]), 
    resize(range(mPromoterGR[mPromoterGR$symbol == x]), width=2000))
  lims <- grLims(gr)
  message(sprintf("Trying %s...", x))
  gr <- GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start=lims[1], end=lims[2])
    )
  gr
  }) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
  width=width(plotRegions) + 0.1*width(plotRegions), 
  fix="center")

plotRegions <- resize(plotRegions, width=ifelse(width(plotRegions) > 100000, width(plotRegions), 100000), fix="center")

fclustOrder <- fineOrder %>% unname() %>% unlist()

p <- plotBrowserTrack(
    ArchRProj = atac_proj, 
    groupBy = "FineNamedClust",
    useGroups = fclustOrder,
    features = flist,
    pal = FineNamedClustCmap,
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(8, 0.8, 1.25, 0.5),
    region = plotRegions, 
    loops = getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff), # All peak-to-gene loops
    tileSize=250,
    minCells=100,
    title = markerGenes
)

plotPDF(p, 
    name = "selected_snp_tracks_100k_width_allClust_allLoops.pdf", 
    ArchRProj = atac_proj, 
    addDOC = FALSE, 
    width = 8, height = 12)

##################################################
# Violin plots of RNA expression for select genes
##################################################

markerGenes <- c("IL33", "SOX8", "HP", "TBX3", "ID2", "FOXN3", "FOXD1", "FOXG1", "FOXl1", "FOXJ2", "FOXO6", "FOXP2", "FOXO1", "NFIA", "NFIC", "FLI1", "ZNF496")

data_mat <- assays(GEmat)[[1]]
rownames(data_mat) <- rowData(GEmat)$name
markerGenes <- markerGenes[which(markerGenes %in% rownames(data_mat))]
sub_mat <- data_mat[markerGenes,]

grouping_data <- data.frame(cluster=factor(atac_proj$FineNamedClust, 
  ordered=TRUE, levels=fclustOrder))
rownames(grouping_data) <- getCellNames(atac_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

pList <- list()
for(gn in markerGenes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  set.seed(1)
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% fclustOrder,] %>% as.data.frame()

  covarLabel <- "cluster"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
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
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Selected_Expression_Violin_gkmsvm_model_genes_FineNamedClust_unclipped.pdf"), width=10, height=4)
pList
dev.off()

