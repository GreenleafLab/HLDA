#!/usr/bin/env Rscript

##################################################
# Create filtered cell barcodes for velocyto
##################################################

suppressPackageStartupMessages({
	library(dplyr)
	library(tidyr)
	library(readr)
	library(Seurat)
  library(reticulate)
	library(sceasy)
})

#Set/Create Working Directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/3_trajectories/adata"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

##### Read in data #####
message("Reading in data...")
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# Add post conception weeks as integer for WOT time series analysis
obj$PCW <- stringr::str_sub(obj$age, start = -2, end = -1) %>% as.integer()

##### Convert to AnnData #####
sceasy::convertFormat(obj, from="seurat", to="anndata", outFile='lda.h5ad')

##### Convert to AnnData with sample name removed from cell barcode and a numerical index added #####
#Convert cell barcode to remove sample info
velocity_obj <- obj
velocity_obj$barcode <- colnames(velocity_obj)
stringr::str_sub(velocity_obj$barcode, 0, -19) <- ""

samples <- unique(velocity_obj$Sample)
suffix <- paste0("-", seq(0,length(samples)-1))

for (i in seq_along(samples)) {
  sample = samples[i]
  end = suffix[i]
  message(paste0("Sample is ", sample))
  message(paste0("Suffix is ", end))
  velocity_obj$barcode[which(velocity_obj$Sample == sample)] <- paste0(velocity_obj$barcode[which(velocity_obj$Sample == sample)], end)
}

sceasy::convertFormat(velocity_obj, from="seurat", to="anndata", outFile='lda_cleaned.h5ad', assay = "RNA")

##### Convert all subclustered objects to AnnData #####
subdir <- "/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/subclustering_final"

clustergroups <- c("Epithelial", "Endothelial", "Immune_v2", "Stromal")

for (cg in clustergroups) {
  so <- readRDS(paste0(subdir, sprintf("/%s/%s.rds", cg, cg)))
  
  # Add PCW info
  so$PCW <- stringr::str_sub(so$age, start = -2, end = -1) %>% as.integer()
  
  #Convert cell barcode to remove sample info
  velocity_obj <- so
  velocity_obj$barcode <- colnames(velocity_obj)
  stringr::str_sub(velocity_obj$barcode, 0, -19) <- ""
  
  samples <- unique(velocity_obj$Sample)
  suffix <- paste0("-", seq(0,length(samples)-1))
  
  for (i in seq_along(samples)) {
    sample = samples[i]
    end = suffix[i]
    message(paste0("Sample is ", sample))
    message(paste0("Suffix is ", end))
    velocity_obj$barcode[which(velocity_obj$Sample == sample)] <- paste0(velocity_obj$barcode[which(velocity_obj$Sample == sample)], end)
  }
  
  sceasy::convertFormat(velocity_obj, from="seurat", to="anndata", outFile=paste0(cg,'_cleaned.h5ad'), assay = "RNA")
}

#####################################################