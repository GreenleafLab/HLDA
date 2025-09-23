#!/usr/bin/env Rscript

##########################################################################
# Calculate AUROC and AUPRC for cross-validation gkmSVM
##########################################################################

# See:
# https://github.com/Dongwon-Lee/lsgkm/
# https://github.com/kundajelab/lsgkm-svr
# https://github.com/kundajelab/gkmexplain/blob/master/dsQTL/gm12878_sequence_sets/compute_auroc.py

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(PRROC)
  library(ComplexHeatmap)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# Read in cross validation results from gkmSVM train
resultDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM/model_predictions"
plotDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM"

# Use capture group to get only files that predict their own cell type
res_files <- list.files(
  path=resultDir, 
  pattern="^(.*\\.)fold.*\\.pred\\.\\1.*\\.txt$", 
  full.names=TRUE
  )

# Files have the format <path>/<cluster>.<fold#>.pred.<target_cell_type>.txt
clusters <- str_replace(basename(res_files), "\\.fold[0-9]+\\..*\\.txt$", "")
fold <- str_extract(basename(res_files), "\\.fold[0-9]+\\.") %>% str_extract(.,"[0-9]+") %>% as.numeric()
truenull <- str_match(basename(res_files), "\\.(true|null)\\.txt$")[,2]
meta_df <- data.frame(cluster=clusters, fold=fold, type=truenull, res_file=res_files)


# Calculate AUROC and AUPRC for each prediction:
getPerformancefromFiles <- function(types, files){
  # Read in gkmSVM cv files and calculate AUROC
  file_df <- data.frame(type=types, res_file=files)
  true_dt <- fread(file_df[file_df$type == "true","res_file"])
  null_dt <- fread(file_df[file_df$type == "null","res_file"])
  # If null is longer than true, sample null sequences for calculating AUROC
  if(nrow(true_dt) < nrow(null_dt)){
    set.seed(1)
    null_dt <- null_dt[sample(1:nrow(null_dt), size=nrow(true_dt), replace=FALSE)]
  }
  AUROC <- PRROC::roc.curve(scores.class0 = true_dt[[2]], scores.class1 = null_dt[[2]])$auc
  AUPRC <- PRROC::pr.curve(scores.class0 = true_dt[[2]], scores.class1 = null_dt[[2]], dg.compute=FALSE)$auc.integral
  data.frame("AUROC"=AUROC, "AUPRC"=AUPRC)
}

summary_df <- meta_df %>% group_by(cluster, fold) %>% do(getPerformancefromFiles(as.character(.$type), as.character(.$res_file))) %>% as.data.frame()

# Specify order of clusters (FineNamedClust)
clustOrder <- fineOrder %>% unname() %>% unlist()

# Only keep cluster names that were sufficient size for training
clustOrder <- clustOrder[clustOrder %in% summary_df$cluster]

# Reorder the cluster names
summary_df$cluster <- factor(summary_df$cluster, levels=clustOrder, ordered=TRUE)

# Load color map for FineNamedClust
cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds")

# Plot a violin / box plot
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

# Plot AUROC
p <- (
  ggplot(summary_df, aes(x=cluster, y=AUROC, color=cluster, fill=cluster))
  + geom_boxplot(alpha=0.5, outlier.shape=NA) # Hide fliers (we show them with geom_jitter)
  + geom_jitter(aes(group=cluster), size=0.75, color="black",
     position=position_jitterdodge(seed=1, jitter.width=7.0, jitter.height=0.0, dodge.width=dodge_width))
  + scale_y_continuous(limits=c(0.8,1.0), expand=c(0,0))
  + scale_color_manual(values=cmap)
  + scale_fill_manual(values=cmap)
  + xlab("")
  + ylab("AUROC")
  + theme_BOR(border=TRUE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          legend.position="none", # Remove legend
          axis.text.x = element_text(angle=90, hjust=1)) 
)

pdf(paste0(plotDir, "/cv_AUROC_boxplot_1000bpNC_fullLim.pdf"), width=8, height=4)
p
dev.off()

# Plot AUPRC
p <- (
  ggplot(summary_df, aes(x=cluster, y=AUPRC, color=cluster, fill=cluster))
  + geom_boxplot(alpha=0.5, outlier.shape=NA) # Hide fliers (we show them with geom_jitter)
  + geom_jitter(aes(group=cluster), size=0.75, color="black",
     position=position_jitterdodge(seed=1, jitter.width=7.0, jitter.height=0.0, dodge.width=dodge_width))
  + scale_y_continuous(limits=c(0.8,1.0), expand=c(0,0))
  + scale_color_manual(values = cmap)
  + scale_fill_manual(values = cmap)
  + xlab("")
  + ylab("AUPRC")
  + theme_BOR(border=TRUE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          legend.position="none", # Remove legend
          axis.text.x = element_text(angle=90, hjust=1)) 
)

pdf(paste0(plotDir, "/cv_AUPRC_boxplot_1000bpNC_fullLim.pdf"), width=8, height=4)
p
dev.off()

###################################################################################