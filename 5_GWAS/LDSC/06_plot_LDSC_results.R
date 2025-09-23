#!/usr/bin/env Rscript

##########################################################################
# Make plots of LDSC results
##########################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(reshape2)
  library(ComplexHeatmap)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# Load colormaps
broadClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/lungClusterColors.rds") %>% unlist()
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds") %>% unlist()

# Identify all result files
resultDir <- "/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/LDSR/FineNamedClust_h2_results"
setwd(resultDir)
resultFiles <- list.files(path=resultDir, pattern=".results$", full.names=TRUE)

readH2results <- function(filename){
  # Read in h2 *.results file
  splits <- strsplit(basename(filename), split="\\.")[[1]]
  ct <- splits[1]
  gwas <- paste0(splits[2:(length(splits)-1)], collapse=".")
  # The first row corresponds to the category of interest, with the remaining being components
  # of the baseline model
  res <- fread(filename)[1,] %>% unlist()
  res[1] <- ct
  res <- c(res, "gwas"=gwas)
  res
}

# Do not use 'allPeaks' or clusters that had <10000 specific peaks
exclude <- c("allPeaks")

results <- lapply(resultFiles, readH2results) %>% do.call(rbind,.) %>% as.data.frame()
results <- results[results$Category %ni% exclude,]
results[,2:7] <- sapply(results[,2:7], as.numeric)
results$Enrichment_FDR <- p.adjust(results$Enrichment_p, method="fdr")
results <- results[order(results$Enrichment_FDR, decreasing=FALSE),]

# Heatmap of enrichment with *** indicating FDR thresholds
# Specify order of clusters (Fine Clust)
clustOrder <- fineOrder %>% unname() %>% unlist()
clustOrder <- clustOrder[clustOrder %in% unique(results$Category)]

mat <- unmelt(results, row_col="gwas", col_col="Category", val_col="Enrichment")
pmat <- unmelt(results, row_col="gwas", col_col="Category", val_col="Enrichment_FDR")

# Separate negative traits (messes with clustering)
neg_mat <- mat[c(
  "PASS_FetalBirthWeight_Warrington2019",
  "PASS_Autism", 
  "PASS_Schizophrenia", 
  "PASS_Neuroticism", 
  "PASS_Anorexia", 
  "PASS_BMI1", 
  "PASS_Years_of_Education1", 
  "UKB_460K.blood_RED_COUNT"
  ), clustOrder]
pos_mat <- mat[rownames(mat) %ni% rownames(neg_mat),]

pos_mat <- pos_mat[,clustOrder]
pos_mat <- prettyOrderMat(pos_mat, clusterCols=FALSE, cutOff=1.0)$mat
mat <- rbind(pos_mat, neg_mat)
pmat <- pmat[rownames(mat), colnames(mat)]

# Get colors for each cluster annotation
#inv.FineClust <- invertList(atac.FineClust)
#broadClust <- unlist(inv.FineClust)[clustOrder] %>% gsub('[0-9]+', '', .) %>% sub('.', '', .)
colors <- FineNamedClustCmap[clustOrder]
names(colors) <- clustOrder

pdf("/oak/stanford/groups/wjg/skim/projects/LDA/GWAS/LDSR/FineNamedClust_h2_results_hm_FineClust.pdf", width=15, height=8)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=colnames(mat),col=list(atac_cluster=colors), 
  show_legend=c("atac_cluster"=FALSE))
hm <- BORHeatmap(
  mat, 
  limits=c(-75,75), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$brewer_yes,
  top_annotation = ta,
  row_names_side = "left",
  width = ncol(mat)*unit(0.5, "cm"),
  height = nrow(mat)*unit(0.5, "cm"),
  legendTitle="Enrichment",
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  cell_fun=function(j,i,x,y,w,h,col){
    if(pmat[i,j] < 0.0005){
      grid.text("***",x,y)
    }else if(pmat[i,j] < 0.005){
      grid.text("**",x,y)
    }else if(pmat[i,j] < 0.05){
      grid.text("*",x,y)
    }else{
      grid.text("",x,y)
    }
  }
  )
draw(hm)
dev.off()

##########################################################################



clustOrder <- c(
    #Epithelial
    "APr1", #
    "APr2", #
    "APr3", #
    "APr4", #
    "PNEC1", # ASCL1+ NEUROD1+
    "PNEC2", # GRP+ PNEC
    "PNEC3", # GHRL+ PNEC
    "eCili", #
    "lCili", #
    "EpiC", #EpiC = cycling epithelial cells
    "TiP1", #
    "TiP2", #
    "AT2l", #
    "AT1l1", # 
    "AT1l2", #
    "AT1l3", 
    "Meso", #
    "Schw", # CDH19+ Schwann
    #Stromal
    "eAlvF", #
    "lAlvF", #
    "MyoF", #
    "AdvF", #
    "aSMC", #
    "ePeri", #SULT1E1+ pericytes (ePeri)
    "lPeri", #LRRTM4+, RBFOX1+, PAG1+, pericytes (lPeri)
    "vSMC1", 
    "vSMC2", #CSMD1+, PLN+, RCAN2+ vSMC (vSMC2)
    "Chdr", #
    #Endothelial
    "egCap", # early gCap
    "lgCap", # Late gCap
    "eAero",
    "lAero", # Aero = Aerocytes
    "eArtr", #
    "lArtr", #arteries SERPINE2+, FBLN5+
    "Veno", #
    "Lymp", #
    #Immune
    "Mono", # (Monocyte VCAN+, AQP9+), 
    "Dc", #(Dendritic cell HLA-DPB1, DPA1 DRB1 HLA positive)
    "Bc", #(Bcells JCHAIN+, IGHM+ BANK1+ )
    "Neut", # Neutrophils LTF+, DEFA3+ (secondary granules of neutrophils and cytotoxic peptides in neutrophil granules) 
    "NK", #
    "Tc" #
    )
