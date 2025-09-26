#!/usr/bin/env Rscript

##########################################################################
# Analysis using finemapped SNPs
##########################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(plyranges)
  library(data.table)
  library(stringr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(parallel)
  library(ggrepel)
  library(ComplexHeatmap)
})

# Set Threads to be used
ncores <- 12

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# set working directory (The directory of the full preprocessed archr project)
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda"
fm_dir <- "/oak/stanford/groups/wjg/skim/projects/LDA/resources/gwas/PICS2"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

topN <- 80 # Number of genes to plot in heatmap

##########################################################################################
# Preparing Data
##########################################################################################

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")
plotDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/finemapping_analysis"
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

raw_finemapped_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))

# Some of the fine-mapped SNPs are duplicated (i.e. the Finacune SNPs sometimes have both FINEMAP and SuSiE finemapping results)
# Deduplicate trait-SNP pairs prior to proceeding with enrichment analyses:
raw_finemapped_gr <- raw_finemapped_gr[order(raw_finemapped_gr$fm_prob, decreasing=TRUE)]
raw_finemapped_gr$trait_snp <- paste0(raw_finemapped_gr$disease_trait, "_", raw_finemapped_gr$linked_SNP)
raw_finemapped_gr <- raw_finemapped_gr[!duplicated(raw_finemapped_gr$trait_snp)] %>% sort()

allGenesGR <- getGenes(atac_proj)

# P2G definition cutoffs
corrCutoff <- 0.45       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Exporting disease traits by source and frequency observed for selecting traits
disease_traits <- raw_finemapped_gr$disease_trait %>% table() %>% as.data.frame()

colnames(disease_traits) <- c("disease_trait", "Freq")
source_traits <- data.frame(disease_trait = raw_finemapped_gr$disease_trait, source = raw_finemapped_gr$source) %>% unique()

disease_traits <- left_join(disease_traits, source_traits, by = "disease_trait")
readr::write_csv(disease_traits, "/oak/stanford/groups/wjg/skim/projects/LDA/resources/gwas/disease_traits.csv")

##########################################################################################
# Plot finemapping probability by peak overlaps
##########################################################################################

prob_bins <- c(0, 0.01, 0.05, 0.1, 0.25, 1)

getPeakOLpct <- function(full_gr, peaks_gr, breaks=c(0, 0.05, 0.1, 0.25, 1)){
  # Calculate percent of SNPs overlapping peaks by bins of finemapped probability
  bin_ids <- makeBins(full_gr$fm_prob, breaks=breaks)$ids
  peak_ol_pct <- sapply(1:(length(breaks)-1), function(b){
    gr <- full_gr[bin_ids == b]
    nol <- overlapsAny(gr, peaks_gr) %>% sum()
    nol/length(gr)
  })
  bin_names <- sapply(2:length(breaks), function(i) paste0(breaks[i-1], "-", breaks[i]))
  names(peak_ol_pct) <- bin_names
  # Get number of SNPs in each bin
  freqs <- getFreqs(bin_ids)
  freqs <- freqs[order(as.integer(names(freqs)))]
  names(freqs) <- bin_names
  list(ol_pct=peak_ol_pct, ol_nSNPs=freqs)
}

disease_traits <- list(
  lung_diseases_a=c(
    "Asthma",
    "Asthma (childhood onset)",
    "Asthma (adult onset)",
    "Asthma onset (childhood vs adult)",
    "Interstitial lung disease",
    "Lung disease severity in cystic fibrosis",
    "CFTR mutation F508del heterozygosity in cystic fibrosis",
    "Chronic obstructive pulmonary disease",
    "Pulmonary fibrosis",
    "Bronchopulmonary dysplasia",
    "Pneumonia",
    "Atopic asthma"
    ),
  lung_diseases_b=c(
    "Childhood asthma with severe exacerbations",
    "Asthma (age of onset)",
    "Asthma (time to onset)",
    "Asthma (moderate or severe)",
    "Cystic fibrosis severity",
    "Chronic obstructive pulmonary disease in non-current smokers",
    "Early chronic obstructive pulmonary disease in never smokers",
    "Idiopathic pulmonary fibrosis",
    "Sarcoidosis (non-Lofgren's syndrome without extrapulmonary manifestations)",
    "COVID-19 (hospitalized vs population)",
    "EGFR mutation-positive lung adenocarcinoma",
    "Lung cancer",
    "Post bronchodilator FEV1/FVC ratio",
    "Chronic obstructive pulmonary disease liability (machine learning-based score)"
    ),
  lung_cancers=c(
    "Lung adenocarcinoma",
    "Squamous cell lung carcinoma",
    "Familial squamous cell lung carcinoma",
    "Non-small cell lung cancer",
    "Small cell lung carcinoma",
    "Familial lung adenocarcinoma"
    ),
  lung_phenotypes=c(
    "FEV1FVC",
    "FEV1",
    "FVC",
    "Lung function (FEV1)",
    "Lung function (FEV1/FVC)",
    "Lung function (FVC)",
    "Peak expiratory flow",
    "Pulmonary function", 
    "Post bronchodilator FEV1 in COPD"
    ),
  brain_traits=c(
    "Major depressive disorder",
    "Schizophrenia",
    "Parkinson's disease",
    "Neuroticism"
  ),
  other_traits=c(
    "BMI", # Finacune finemapped
    "Educational attainment (years of education)",
    "Red blood cell count",
    "Height",
    "Balding_Type4",
    "SBP" # Finacune finemapped
  )
)
all_disease_traits <- do.call(c,disease_traits) %>% unname()

pList <- lapply(all_disease_traits, function(dt){
  message(sprintf("Testing trait %s...", dt))
  ol <- getPeakOLpct(raw_finemapped_gr[raw_finemapped_gr$disease_trait == dt], peaks_gr=allPeaksGR, breaks=prob_bins)
  pct_ol <- ol$ol_pct
  nSNPs <- ol$ol_nSNPs
  xlabs <- paste0(names(nSNPs), sprintf("\n(%s)", nSNPs))
  df <- data.frame(fm_probs=names(pct_ol), percent_snps_overlapping=pct_ol)
  df$fm_probs <- factor(df$fm_probs, levels=names(pct_ol), ordered=TRUE)
  qcBarPlot(df, cmap="grey80", barwidth=0.9, border_color="black") + 
    scale_y_continuous(limits=c(0, 0.4), expand = c(0, 0)) + 
    ggtitle(dt) + scale_x_discrete(labels=xlabs)
  })

pdf(paste0(plotDir, "/pctPeaksOL_by_fm_probs_eachTrait.pdf"), width=5, height=5)
pList
dev.off()

disease_traits <- list(
  lung_diseases_a=c(
    "Asthma",
    "Asthma (childhood onset)",
    "Asthma (adult onset)",
    "Asthma onset (childhood vs adult)",
    "Bronchopulmonary dysplasia",
    "Atopic asthma",
    "Asthma (age of onset)",
    "Early chronic obstructive pulmonary disease in never smokers",
    "Chronic obstructive pulmonary disease liability (machine learning-based score)"
    ),
  lung_diseases_b=c(
    "Interstitial lung disease",
    "Lung disease severity in cystic fibrosis",
    "CFTR mutation F508del heterozygosity in cystic fibrosis",
    "Chronic obstructive pulmonary disease",
    "Pulmonary fibrosis",
    "Pneumonia",
    "Childhood asthma with severe exacerbations",
    "Asthma (time to onset)",
    "Asthma (moderate or severe)",
    "Cystic fibrosis severity",
    "Chronic obstructive pulmonary disease in non-current smokers",
    "Idiopathic pulmonary fibrosis",
    "Sarcoidosis (non-Lofgren's syndrome without extrapulmonary manifestations)",
    "COVID-19 (hospitalized vs population)",
    "EGFR mutation-positive lung adenocarcinoma",
    "Lung cancer",
    "Post bronchodilator FEV1/FVC ratio"
    ),
  lung_cancers=c(
    "Lung adenocarcinoma",
    "Squamous cell lung carcinoma",
    "Familial squamous cell lung carcinoma",
    "Non-small cell lung cancer",
    "Small cell lung carcinoma",
    "Familial lung adenocarcinoma"
    ),
  lung_phenotypes_a=c(
    "FEV1FVC",
    "Lung function (FEV1/FVC)",
    "Peak expiratory flow"
    ),
  lung_phenotypes_b=c(
    "FEV1",
    "Lung function (FEV1)",
    "Lung function (FVC)",
    "Pulmonary function", 
    "Post bronchodilator FEV1 in COPD"
    ),
  brain_traits=c(
    "Major depressive disorder",
    "Schizophrenia",
    "Parkinson's disease",
    "Neuroticism"
  ),
  other_traits=c(
    "BMI", # Finacune finemapped
    "Educational attainment (years of education)",
    #"Red blood cell count",
    #"Height",
    #"Balding_Type4",
    "SBP" # Finacune finemapped
  )
)

# Plot box plot showing increase in fraction of snps overlapping peaks in large groups of SNP categories
use_groups <- c("lung_diseases_a", "lung_phenotypes_a", "brain_traits", "other_traits")
sub_traits_list <- disease_traits[use_groups]

ol_df <- lapply(names(sub_traits_list), function(x){
  message(sprintf("Calculating overlaps for group %s...", x))
  dts <- sub_traits_list[[x]]
  df <- lapply(dts, function(dt){
      ol <- getPeakOLpct(raw_finemapped_gr[raw_finemapped_gr$disease_trait == dt], peaks=allPeaksGR, breaks=prob_bins)
      pct_ol <- ol$ol_pct
      nSNPs <- ol$ol_nSNPs
      df <- data.frame(fm_probs=names(pct_ol), nSNPs=nSNPs, pct_snps_ol=pct_ol)
      df$fm_probs <- factor(df$fm_probs, levels=names(pct_ol), ordered=TRUE)
      df$trait <- dt
      df
    }) %>% do.call(rbind,.)
  df$trait_set <- x
  df
  }) %>% do.call(rbind,.)

ol_df$trait_set <- factor(ol_df$trait_set, ordered=TRUE, levels=use_groups)

# Plot multiple category box plot
cmap <- cmaps_BOR$stallion
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

p <- (
  ggplot(ol_df, aes(x=fm_probs, y=pct_snps_ol, color=trait_set, fill=trait_set))
  + geom_boxplot(alpha=0.5, outlier.shape = NA) # Hide fliers (we show them with geom_jitter)
  + geom_jitter(aes(group=trait_set), size=1.0, color="black",
     position=position_jitterdodge(seed=1, jitter.width=0.25, jitter.height=0.0, dodge.width=dodge_width))
  + scale_fill_manual(values=cmap)
  + scale_color_manual(values=cmap)
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
          axis.text.x = element_text(angle = 90, hjust = 1)) 
  + scale_y_continuous(limits=c(0, 0.39), expand = c(0, 0))
)
pdf(paste0(plotDir, "/pctPeaksOL_by_fm_probs_og.pdf"), width=10, height=4)
p
dev.off()


# Read in cell type specific peaks
cts_peak_files <- list.files(
  path="/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/LDSR/FineNamedClust_specific_peaks",
  pattern="*.bed$",
  full.names=TRUE
)
names(cts_peak_files) <- basename(cts_peak_files) %>% sapply(., function(x) gsub("_specific_peaks.bed", "", x)) %>% unname()

cts_peaks <- lapply(cts_peak_files, function(pf){
  dt <- fread(pf)
  GRanges(seqnames=dt$V1, IRanges(dt$V2, end=dt$V3))
  })
names(cts_peaks) <- names(cts_peak_files)

fisherTestSNPs <- function(peaks_gr, snps_gr, disease_trait){
  # Calculate fisher enrichment for a specific trait in a given peak set
  ######################################################################
  # peaks_gr = peakset to use for testing enrichment
  # snps_gr = full fine-mapping GR
  # disease_trait = which trait to test enrichment for

  # SNPs overlapping peakset
  nondis_gr <- snps_gr[snps_gr$disease_trait != disease_trait]
  nol <- overlapsAny(nondis_gr, peaks_gr) %>% sum() # n non-disease SNPs overlapping
  nnol <- length(nondis_gr) - nol # n non-disease SNPs not overlapping

  # Trait SNPs overlapping peakset
  dis_gr <- snps_gr[snps_gr$disease_trait == disease_trait] 
  dnol <- overlapsAny(dis_gr, peaks_gr) %>% sum() # n disease SNPs overlapping
  dnnol <- length(dis_gr) - dnol # n disease SNPs not overlapping

  OR <- (dnol/dnnol)/(nol/nnol)
  pval <- fisher.test(matrix(c(dnol, dnnol, nol, nnol),2,2), alternative="greater")$p.value

  list(trait=disease_trait, ol_dis_snps=dnol, nol_dis_snps=dnnol, ol_snps=nol, nol_snps=nnol, OR=OR, fisher_pval=pval)
}

# (Both PICS and finacune finemapped SNPs have ~20% with prob > 0.05) Previously tried 0.01
filt_snp_gr <- raw_finemapped_gr[raw_finemapped_gr$fm_prob >= 0.05]

res_df <- lapply(names(cts_peaks), function(ct){
    message(sprintf("Testing %s peaks...", ct))
    peaks_gr <- cts_peaks[[ct]]
    dis_res <- lapply(all_disease_traits, function(dt){
      fisherTestSNPs(peaks_gr, filt_snp_gr, dt)
    })
    dis_res <- do.call(rbind, lapply(dis_res, data.frame))
    dis_res$cluster <- ct
    dis_res
  }) %>% do.call(rbind,.)

res_df$padj <- p.adjust(res_df$fisher_pval, method="fdr")
res_df <- res_df[order(res_df$padj, decreasing=FALSE),]
res_df$mlog10padj <- -log10(res_df$padj)

saveRDS(res_df, paste0(fm_dir, "/res_df_SNPEnrichment.rds"))

# Plot Dot Plot of fischer enrichment
# Specify order of clusters (CustomNamedClust)
clustOrder <- unname(unlist(fineOrder))

clustOrder <- clustOrder[clustOrder %in% unique(res_df$cluster)]

# Plot separate dot plot for each group:
max_mlogpval <- 15
colorLims <- c(0, max_mlogpval)
sizeLims <- c(min(res_df$OR), max(res_df$OR))

pList <- list()
for(grp in names(disease_traits)){

  sub_res_df <- res_df[res_df$trait %in% disease_traits[[grp]],]
  # Determine cluster and gene order:
  wide_df <- unmelt(sub_res_df, row_col="trait", col_col="cluster", val_col="OR")
  wide_df <- prettyOrderMat(wide_df[,clustOrder], clusterCols=FALSE)$mat

  # Prepare for plotting
  plot_df <- sub_res_df[sub_res_df$cluster %in% clustOrder,]
  plot_df$grp <- plot_df$cluster
  plot_df$mlog10padj <- ifelse(plot_df$mlog10padj > max_mlogpval, max_mlogpval, plot_df$mlog10padj)

  grp_order <- colnames(wide_df)
  trait_order <- rownames(wide_df) %>% rev()

  pList[[grp]] <- dotPlot(plot_df, xcol="grp", ycol="trait", color_col="mlog10padj", size_col="OR", 
    xorder=unlist(clustOrder), yorder=trait_order, 
    cmap=cmaps_BOR$wolfgang_extra, aspectRatio=nrow(wide_df)/ncol(wide_df), 
    sizeLims=sizeLims, colorLims=colorLims)
}

pdf(paste0(plotDir, "/fmSNP_enrichment_fisher_clusters_dotPlots_og.pdf"), width=15, height=8)
pList
dev.off()

rm(raw_finemapped_gr); gc() 


