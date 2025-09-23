#!/usr/bin/env Rscript

#####################################################################
# Call peaks on clusters from overall clustering
#####################################################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# Working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess"
setwd(wd)

#Load Genome Annotations
addArchRGenome("hg38")
pointSize <- 0.5

##########################################################################################
# Load Previously Prepared ArchR project
##########################################################################################
#proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/Figure01_Data_Summary/integrated_lda")
#saveArchRProject(proj, outputDirectory = paste0(wd, "/lda_v2"),  load = FALSE)
proj <- loadArchRProject(paste0(wd, "/lda_v2"), force = T)

##########################################################################################
# Subcluster each compartment
##########################################################################################

compartments <- list(
    # "Epithelial" = c("APr1", "APr2", "APr3", "APr4", 
    #                  "PNEC1", "PNEC2", "PNEC3", 
    #                  "eCili", "lCili", "TiP1", "TiP2", "AT2l", "AT1l", "EpiC"),
    # "Endothelial" = c("egCap", "lgCap", "eAero", "lAero", "eArtr", "lArtr", "Veno", "Lymp"),
    "Stromal" = c("eAlvF", "lAlvF", "AdvF", "MyoF", "aSMC", "ePeri", "lPeri", "vSMC1", "vSMC2", "Chdr", "Meso", "Schw")#,
    # "Immune" = c("Tc", "NK", "Bc", "Plasma", "MPP", "Neut", "Dc", "Mono1", "Mono2", "IM"),
    # "Alveolar" = c("APr1", "APr2", "APr3", "APr4", "TiP1", "TiP2", "AT2l", "AT1l"),
    # "PNECs" = c("APr2", "APr3", "PNEC1", "PNEC2", "PNEC3"),
    # "Vascular" = c("egCap", "lgCap", "eAero", "lAero", "eArtr", "lArtr", "Veno"),
    # "Stromal_a" = c("eAlvF", "lAlvF", "AdvF", "MyoF", "aSMC"),
    # "Stromal_b" = c("ePeri", "lPeri", "vSMC1", "vSMC2"),
    # "Myeloid" = c("MPP", "Neut", "Dc", "Mono1", "Mono2", "IM"),
    # "Lymphoid" = c("Tc", "NK", "Bc", "Plasma")
    )

subClusterCells <- lapply(compartments, function(x){
    getCellNames(proj)[proj@cellColData$FineNamedClust %in% x]
})

subClusterArchR <- function(proj, subCells, outdir){
  # Subset an ArchR project for focused analysis

  message(sprintf("Subgroup has %s cells.", length(subCells)))
  sub_proj <- subsetArchRProject(
      ArchRProj = proj,
      cells = subCells,
      outputDirectory = outdir,
      dropCells = TRUE,
      force = TRUE
  )
  saveArchRProject(sub_proj)

  # Delete stuff we don't want to copy...
  unlink(paste0(outdir, "/Plots/*"), recursive = TRUE)
  unlink(paste0(outdir, "/Peak2GeneLinks/*"), recursive = TRUE)
}

# Generate subprojects from each compartment
sg_list <- names(compartments)

dir.create("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering", showWarnings = FALSE, recursive = TRUE)

sub_proj_list <- lapply(sg_list, function(sg){
  message(sprintf("Subsetting %s...", sg))
  outdir <- sprintf("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/%s", sg)
  subClusterArchR(proj, subCells=subClusterCells[[sg]], outdir=outdir)
})

