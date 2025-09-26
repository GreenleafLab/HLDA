#!/usr/bin/env Rscript

##################################################
# Create cell sets for WOT 
##################################################

suppressPackageStartupMessages({
	library(dplyr)
	library(tidyr)
	library(readr)
	library(Seurat)
})

#Set/Create Working Directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/5_wot"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

##### Read in data #####
message("Reading in data...")
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1b_rna_preprocess/lda.rds")

##### Function for writing GMT files #####

write_gmt <- function(lst, filename) {
  # assumes that each element of the list will have the fields
  # head, desc, entry
  if (file.exists(filename)) {
    message(paste(filename, "exists, deleting..."))
    file.remove(filename)
  }

  for (i in seq_along(lst)) {
    el <- lst[[i]]
    ncolumns <- 2 + length(el$entry)
    write(c(el$head, el$desc, el$entry), file=filename, sep="\t",
          append=TRUE, ncolumns=ncolumns)
  }
}

##### Create a GMT file for cell clusters #####
# GMT file for cell clusters has the structure below:
# First column: Cell Cluster name
# Second column: Description but can be left blank or NA
# Third column onwards: Cell barcodes in clusters

# This file format is required for WOT

lst <- list()
ClusterNames <- unique(obj$CustomNamedClust)
Idents(obj) <- "CustomNamedClust"

# For each cluster, grab cell barcode and append to outlist
for (ClusterName in ClusterNames) {
  CBs <- WhichCells(obj, idents = ClusterName)
  list <- list(head = ClusterName, desc = NA, entry = CBs) # NA is because desc field can be empty
  lst <- append(lst, list(list)) #append to a list of lists
}

# Save file 
write_gmt(lst, paste0(wd, "/cell_sets_lda.gmt"))

#####################################################