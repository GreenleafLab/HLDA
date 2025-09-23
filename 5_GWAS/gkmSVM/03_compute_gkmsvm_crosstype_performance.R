#!/usr/bin/env Rscript

##########################################################################
# Calculate AUROC for cross-validation gkmSVM
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
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# Read in cross validation results from gkmSVM train
resultDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM/model_predictions"
plotDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/gkmSVM"

res_files <- list.files(
  path=resultDir, 
  pattern="\\.fold0\\.pred\\..*\\.txt$", 
  full.names=TRUE
  )

# Files have the format <path>/<pred_cluster>.pred.<target_type>.[true/null].txt
clusters <- str_replace(basename(res_files), "\\.fold[0-9]+\\..*\\.txt$", "")
targets <- str_match(basename(res_files), "\\.pred\\.([^\\.]*)\\..*")[,2]
truenull <- str_match(basename(res_files), "\\.(true|null)\\.txt$")[,2]

meta_df <- data.frame(cluster=clusters, target=targets, type=truenull, res_file=res_files)

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

# Fairly slow... (15min)
summary_df <- meta_df %>% group_by(cluster, target) %>% do(getPerformancefromFiles(.$type, .$res_file)) %>% as.data.frame()

auroc_mat <- unmelt(summary_df, row_col="cluster", col_col="target", val_col="AUROC")
auprc_mat <- unmelt(summary_df, row_col="cluster", col_col="target", val_col="AUPRC")
# mat <- prettyOrderMat(mat,clusterCols=TRUE, cutOff=1.0)$mat

# Specify order of clusters (fineOrder)
clustOrder <- fineOrder %>% unname() %>% unlist()

# Only keep cluster names that were sufficient size for training
clustOrder <- clustOrder[clustOrder %in% summary_df$cluster]

# AUROC heatmap
pdf(paste0(plotDir, "/crosstype_AUROC_hm_1000bpNC.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  auroc_mat[clustOrder, clustOrder], 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  column_title="Target Cell Type", 
  column_title_side="bottom",
  row_title="Model Cell Type",
  dataColors = cmaps_BOR$solar_extra,
  row_names_side = "left",
  width = ncol(auroc_mat)*unit(0.5, "cm"),
  height = nrow(auroc_mat)*unit(0.5, "cm"),
  legendTitle="AUROC",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()

# AUPRC heatmap
pdf(paste0(plotDir, "/crosstype_AUPRC_1000bpNC.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  auprc_mat[clustOrder, clustOrder], 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  column_title="Target Cell Type", 
  column_title_side="bottom",
  row_title="Model Cell Type",
  dataColors = cmaps_BOR$solar_extra,
  row_names_side = "left",
  width = ncol(auprc_mat)*unit(0.5, "cm"),
  height = nrow(auprc_mat)*unit(0.5, "cm"),
  legendTitle="AUPRC",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
draw(hm)
dev.off()




