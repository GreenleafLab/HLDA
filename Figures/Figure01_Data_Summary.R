#!/usr/bin/env Rscript

########################################
# Prepare figure panels for figure 1
########################################

#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(BuenColors)
})

# Get additional functions, etc.:
scriptPath <- "/oak/stanford/groups/wjg/skim/projects/LDA/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

# set working directory
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/Figure01_Data_Summary"

#Set/Create Working Directory to Folder
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

# set plot directory
plotDir <- paste0(wd)

# Misc options
addArchRGenome("hg38")
pointSize <- 0.3

##########################################################################################
# Preparing Data
##########################################################################################
# Load data objects
atac_proj <- loadArchRProject("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/lda_v2")
rna_proj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")

# Color Maps
compartmentCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/lungClusterColors.rds") %>% unlist()
FineNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_FineNamedClust_cmap.rds") %>% unlist()
BroadNamedClustCmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/scRNA_BroadNamedClust_cmap.rds") %>% unlist()

sample_cmap <- readRDS(paste0("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/sample_cmap.rds"))
gest_age_cmap <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/cmaps/gest_age_cmap.rds")

##########################################################################################
# UMAP of ATAC and RNA with BroadClust and FineClust
##########################################################################################
umapPlots <- list()

plot.name <- "UMAP_ATAC_RNA_FineBroadNamedClust"
dir.name <- paste0(wd, "/", plot.name)
dir.create(dir.name, showWarnings = FALSE, recursive = TRUE)

# ATAC BroadNamedClust on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="BroadNamedClust", embeddingName="UMAP_ATAC")
readr::write_tsv(umapDF, file = paste0(dir.name, "/BroadNamedClust_on_ATAC.tsv"))
umapPlots[["BroadNamedClust_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=BroadNamedClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="BroadNamedClust_on_ATAC", useRaster=TRUE)

# ATAC FineNamedClust on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="FineNamedClust", embeddingName="UMAP_ATAC")
readr::write_tsv(umapDF, file = paste0(dir.name, "/FineNamedClust_on_ATAC.tsv"))
umapPlots[["FineNamedClust_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=FineNamedClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="FineNamedClust_on_ATAC", useRaster=TRUE)

# RNA BroadNamedClust on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$BroadNamedClust)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
readr::write_tsv(umapDF, file = paste0(dir.name, "/BroadNamedClust_on_RNA.tsv"))
umapPlots[["BroadNamedClust_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=BroadNamedClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="BroadNamedClust_on_RNA", useRaster=TRUE)

# RNA FineNamedClust on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$FineNamedClust)
set.seed(1)
readr::write_tsv(umapDF, file = paste0(dir.name, "/FineNamedClust_on_RNA.tsv"))
umapDF <- umapDF[sample(nrow(umapDF)),]
umapPlots[["FineNamedClust_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=FineNamedClustCmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="FineNamedClust_on_RNA", useRaster=TRUE)

pdf(paste0(plotDir,"/",plot.name,".pdf"), width=10, height=8)
umapPlots
dev.off()

##########################################################################################
# UMAP of ATAC and RNA with gestational age and sample identity
##########################################################################################
umapPlots <- list()

plot.name <- "clusters_UMAP_ATAC_RNA_Age_Samples"
dir.name <- paste0(wd, "/", plot.name)
dir.create(dir.name, showWarnings = FALSE, recursive = TRUE)

# Age on RNA clustering:
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$age)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
readr::write_tsv(umapDF, file = paste0(dir.name, "/Ages_on_RNA.tsv"))
umapPlots[["Ages_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=gest_age_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="Ages_on_RNA", useRaster=TRUE)

# Sample on RNA clustering
umapDF <- data.frame(Embeddings(object = rna_proj, reduction = "umap"), rna_proj$Sample)
set.seed(1)
umapDF <- umapDF[sample(nrow(umapDF)),]
readr::write_tsv(umapDF, file = paste0(dir.name, "/Samples_on_RNA.tsv"))
umapPlots[["Samples_on_RNA"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=sample_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="Samples_on_RNA", useRaster=TRUE)

# Age on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="age", embeddingName="UMAP_ATAC")
readr::write_tsv(umapDF, file = paste0(dir.name, "/Ages_on_ATAC.tsv"))
umapPlots[["Ages_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=gest_age_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="Ages_on_ATAC", useRaster=TRUE)

# Sample on ATAC clustering:
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="Sample", embeddingName="UMAP_ATAC")
readr::write_tsv(umapDF, file = paste0(dir.name, "/Samples_on_ATAC.tsv"))
umapPlots[["Samples_on_ATAC"]] <- plotUMAP(umapDF, dataType = "qualitative", cmap=sample_cmap, 
  namedColors=TRUE, point_size=pointSize, covarLabel="Samples_on_ATAC", useRaster=TRUE)

pdf(paste0(plotDir,"/",plot.name,".pdf"), width=10, height=8)
umapPlots
dev.off()


#############################################################################
# Stacked Bar Plots of Sample by Cluster
#############################################################################

rna_proj@meta.data$pcw <- gsub("PCW", "", rna_proj@meta.data$age)

# Order clusters based on increasing average PCW of the cluster
broadOrder <- rna_proj@meta.data %>% 
  group_by(BroadNamedClust) %>%
  summarize(avgPCW = mean(as.integer(pcw))) %>%
  arrange(avgPCW) %>% 
  dplyr::select(BroadNamedClust) %>% as.list() %>% unname() %>% unlist()

# broadOrder <- c(
#     "APr", "PNEC", "Cili", "TiP",  "AT2l", "AT1l", "EpiC", "Meso",
#     "gCap", "Aero", "Artr", "Veno", "Lymp",
#     "AlvF", "AdvF", "MyoF", "aSMC", "Peri", "vSMC", "Chdr",
#     "Tc", "NK", "APCs"
#     )

barWidth <- 1

# BroadNamedClust by Age
clustByAge <- fractionXbyY(rna_proj$BroadNamedClust, rna_proj$age, add_total=TRUE, xname="BroadNamedClust", yname="Gestational Age")
clustByAge$BroadNamedClust <- factor(clustByAge$BroadNamedClust, levels = c(broadOrder, "total"), ordered = TRUE)
pdf(paste0(plotDir, "/barplot_BroadclustByAge.pdf"), height = 5, width = 10)
stackedBarPlot(clustByAge, cmap=gest_age_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# BroadNamedClust by Sample
clustBySamp <- fractionXbyY(rna_proj$BroadNamedClust, rna_proj$Sample, add_total=TRUE, xname="BroadNamedClust", yname="Sample")
clustBySamp$BroadNamedClust <- factor(clustBySamp$BroadNamedClust, levels = c(broadOrder, "total"), ordered = TRUE)
pdf(paste0(plotDir, "/barplot_BroadclustBySample.pdf"), height = 5, width = 10)
stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

barWidth <- 1

# Order clusters based on increasing average PCW of the cluster

# FineNamedClust by Age
clustByAge <- fractionXbyY(rna_proj$FineNamedClust, rna_proj$age, add_total=TRUE, xname="FineNamedClust", yname="Gestational Age")
clustByAge$FineNamedClust <- factor(clustByAge$FineNamedClust, levels = c(fineOrder, "total"), ordered = TRUE)
pdf(paste0(plotDir, "/barplot_FineclustByAge_ordered.pdf"), height = 5, width = 15)
stackedBarPlot(clustByAge, cmap=gest_age_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# FineNamedClust by Sample
clustBySamp <- fractionXbyY(rna_proj$FineNamedClust, rna_proj$Sample, add_total=TRUE, xname="FineNamedClust", yname="Sample")
clustBySamp$FineNamedClust <- factor(clustBySamp$FineNamedClust, levels = c(fineOrder, "total"), ordered = TRUE)
pdf(paste0(plotDir, "/barplot_FineclustBySample.pdf"), height = 5, width = 5)
stackedBarPlot(clustBySamp, cmap=sample_cmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

#############################################################################
# Stacked Bar Plots of Compartment by Age
#############################################################################

barWidth <- 1.0

# Compartment by Age
AgebyCompartment <- fractionXbyY(rna_proj$age, rna_proj$compartment, add_total=TRUE, xname="Gestational Age", yname="Compartment")
pdf(paste0(plotDir, "/barplot_AgeByCompartment.pdf"), height = 3, width = 8)
stackedBarPlot(AgebyCompartment, cmap=compartmentCmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

#############################################################################
# Stacked Bar Plots of Cluster by Age
#############################################################################

barWidth <- 1.0

# BroadNamedClust by Age
AgebyClust <- fractionXbyY(rna_proj$age, rna_proj$BroadNamedClust, add_total=TRUE, xname="Gestational Age", yname="BroadNamedClust")
pdf(paste0(plotDir, "/barplot_AgeByBroadclust.pdf"), height = 5, width = 10)
stackedBarPlot(AgebyClust, cmap=BroadNamedClustCmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

# FineNamedClust by Age
AgebyClust <- fractionXbyY(rna_proj$age, rna_proj$FineNamedClust, add_total=TRUE, xname="Gestational Age", yname="FineNamedClust")
pdf(paste0(plotDir, "/barplot_AgeByFineclust.pdf"), height = 5, width = 10)
stackedBarPlot(AgebyClust, cmap=FineNamedClustCmap, namedColors=TRUE, barwidth=barWidth)
dev.off()

#############################################################################
# Stacked barplot of clusters by age ordered by mean age
#############################################################################
# FineNamedClust by Age
clustByAge <- fractionXbyY(rna_proj$FineNamedClust, rna_proj$age, add_total=TRUE, xname="FineNamedClust", yname="Gestational Age")
clustByAge$FineNamedClust <- factor(clustByAge$FineNamedClust, levels = c(fineOrder, "total"), ordered = TRUE)

# Add post conception weeks as integer to calculate mean age of each cluster
rna_proj$PCW <- stringr::str_sub(rna_proj$age, start = -2, end = -1) %>% as.integer()
mean.ages <- data.frame(FineNamedClust = rna_proj@meta.data$FineNamedClust,
                        PCW = rna_proj@meta.data$PCW) %>% 
  group_by(FineNamedClust) %>% 
  summarize(mean.age = mean(PCW)) %>%
  arrange(mean.age)

# Save mean age for each cluster
readr::write_tsv(rna_proj@meta.data %>% group_by(FineNamedClust) %>% summarise(mean_age = mean(PCW)), file = paste0(wd, "/mean_age_per_FineNamedClust.tsv"))

# For each compartment order clusters based on mean age
organized_clusters <- list()
for (i in 1:length(fineOrder)) {
  clusters <- fineOrder[i] %>% unlist() %>% unname()
  sub.mean.ages <- mean.ages[which(mean.ages$FineNamedClust %in% clusters),] %>% arrange(mean.age)
  organized_clusters <- c(organized_clusters, unlist(unname(sub.mean.ages$FineNamedClust)))
}
organized_clusters <- unlist(unname(organized_clusters))

# Reorder clusters based on mean average age
clustByAge$FineNamedClust <- factor(clustByAge$FineNamedClust, levels = organized_clusters)
clustByAge$`Gestational Age` <- factor(clustByAge$`Gestational Age`, levels = c("PCW23", "PCW21", "PCW19", "PCW16", "PCW15", "PCW14", "PCW12"))

saveRDS(clustByAge, paste0(plotDir, "/barplot_FineClustByAge_OrderedPerCompartment.rds"))

pdf(paste0(plotDir, "/barplot_FineClustByAge_OrderedPerCompartment.pdf"), h = 3, w = 7)
ggplot(clustByAge, aes(fill=`Gestational Age`, y=proportion, x=FineNamedClust)) + 
  geom_bar(position="fill", stat="identity") +
  scale_x_discrete(labels = levels(clustByAge$FineNamedClust)) +  # To display the reordered labels
  scale_fill_manual(values = gest_age_cmap) +
  theme_BOR() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.title = element_blank()
  )
dev.off()

#############################################################################
# Stacked barplot of Proximal vs Distal normalized for cell counts in each sample
#############################################################################
rna_proj@meta.data <-
  rna_proj@meta.data %>% mutate(location = ifelse(grepl("[PD]$", Sample),
                                                  ifelse(grepl("P$", Sample), "Proximal", "Distal"),
                                                  NA))
sub_rna <- rna_proj[,!is.na(rna_proj$location)]

clustByLoc <- normFractionXbyY(sub_rna$FineNamedClust, sub_rna$location, 
                           xname="FineNamedClust", 
                           yname="Location")

sortByLoc <- clustByLoc %>%
  filter(Location == "Proximal") %>% arrange(proportion)

# For each compartment order clusters based on proportion of proximal to distal
organized_clusters <- list()
for (i in 1:length(fineOrder)) {
  clusters <- fineOrder[i] %>% unlist() %>% unname()
  sub.loc.prop <- sortByLoc[which(sortByLoc$FineNamedClust %in% clusters),] %>% arrange(desc(proportion))
  organized_clusters <- c(organized_clusters, as.vector(sub.loc.prop$FineNamedClust))
}
organized_clusters <- unlist(unname(organized_clusters))

# Reorder clusters based on mean average age
clustByLoc$FineNamedClust <- factor(clustByLoc$FineNamedClust, levels = organized_clusters)

# To show number of cells in eeach cluster
# cluster_count <- sub_rna$FineNamedClust %>% table() %>% as.data.frame()
# colnames(cluster_count) <- c("FineNamedClust", "Freq")
# cluster_count$label <- paste0(cluster_count[,1]," (",cluster_count[,2],")")
# ordered_clust <- data.frame(FineNamedClust = levels(clustByLoc$FineNamedClust))
# cluster_count <- left_join(ordered_clust, cluster_count, by = "FineNamedClust")

pdf(paste0(plotDir, "/barplot_FineClustByLocation_OrderedPerCompartment.pdf"), h = 2, w = 9)
ggplot(clustByLoc, aes(fill=Location, y=proportion, x=FineNamedClust)) + 
  geom_bar(position="fill", stat="identity") +
  scale_x_discrete(labels = levels(clustByLoc$FineNamedClust)) +  # To display the reordered labels
  scale_fill_manual(values = cmaps_BOR$stallion[1:2]) +
  geom_hline(yintercept = 0.5) +
  theme_BOR() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title.x = element_blank()
  )
dev.off()

#########################################################
# Stacked bar plot of age by cluster for each compartment 
#########################################################

Idents(rna_proj) <- "FineNamedClust"

mat.list <- list()
for (i in 1:length(fineOrder)) {
  clusters <- fineOrder[i] %>% unlist() %>% unname()
  name <- names(fineOrder[i])
  
  sub_rna <- subset(rna_proj, idents = clusters)
  ageByClust <- fractionXbyY(sub_rna$age, sub_rna$FineNamedClust, 
                             xname="Age", 
                             yname="FineNamedClust", add_total = T)
  ageByClust$FineNamedClust <- factor(ageByClust$FineNamedClust, levels = clusters)
  mat.list[[name]] <- ageByClust
}

pdf(paste0(plotDir, "/barplot_GestationalAgeByFineNamedClust_PerCompartment.pdf"), h = 4, w = 5)
ggplot(mat.list[["Epithelial"]], aes(fill=FineNamedClust, y=proportion, x=Age)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = FineNamedClustCmap[fineOrder$Epithelial]) +
  theme_BOR() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title.x = element_blank()
  )

ggplot(mat.list[["Endothelial"]], aes(fill=FineNamedClust, y=proportion, x=Age)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = FineNamedClustCmap[fineOrder$Endothelial]) +
  theme_BOR() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title.x = element_blank()
  )

ggplot(mat.list[["Stromal"]], aes(fill=FineNamedClust, y=proportion, x=Age)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = FineNamedClustCmap[fineOrder$Stromal]) +
  theme_BOR() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title.x = element_blank()
  )

ggplot(mat.list[["Immune"]], aes(fill=FineNamedClust, y=proportion, x=Age)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = FineNamedClustCmap[fineOrder$Immune]) +
  theme_BOR() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.title.x = element_blank()
  )
dev.off()


#############################################################################
# QC Plots for RNA
#############################################################################

# Violin plots for RNA QC metrics
p1 <- plotGroups(ArchRProj = atac_proj, groupBy = "age", colorBy = "cellColData", name = "log10(Gex_nUMI)",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = gest_age_cmap)

p2 <- plotGroups(ArchRProj = atac_proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(Gex_nUMI)",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = sample_cmap)

p3 <- plotGroups(ArchRProj = atac_proj, groupBy = "age", colorBy = "cellColData", name = "log10(Gex_nGenes)",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = gest_age_cmap)

p4 <- plotGroups(ArchRProj = atac_proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(Gex_nGenes)",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = sample_cmap)

p5 <- plotGroups(ArchRProj = atac_proj, groupBy = "age", colorBy = "cellColData", name = "Gex_MitoRatio",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = gest_age_cmap)

p6 <- plotGroups(ArchRProj = atac_proj, groupBy = "Sample", colorBy = "cellColData", name = "Gex_MitoRatio",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = sample_cmap)

pdf(paste0(wd, "/QC-RNA-Statistics.pdf"), w = 2.5, h = 2.5)
p1
p2
p3
p4
p5
p6
dev.off()

#plotPDF(p1, p2, p3, p4, p5, p6, name = "QC-RNA-Statistics.pdf", ArchRProj = atac_proj, addDOC = FALSE, width = 4, h = 2)


#############################################################################
# QC Plots for ATAC
#############################################################################

# Plots for ATAC QC metrics
p1 <- plotFragmentSizes(ArchRProj = atac_proj, groupBy = "age", pal = gest_age_cmap)

p2 <- plotFragmentSizes(ArchRProj = atac_proj, groupBy = "Sample", pal = sample_cmap)

p3 <- plotTSSEnrichment(ArchRProj = atac_proj, groupBy = "age", pal = gest_age_cmap)

p4 <- plotTSSEnrichment(ArchRProj = atac_proj, groupBy = "Sample", pal = sample_cmap)

p5 <- plotGroups(ArchRProj = atac_proj, groupBy = "age", colorBy = "cellColData", name = "FRIP",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = gest_age_cmap)

p6 <- plotGroups(ArchRProj = atac_proj, groupBy = "Sample", colorBy = "cellColData", name = "FRIP",
                 plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE, pal = sample_cmap)

pdf(paste0(wd, "/QC-ATAC-Statistics_FragTSS.pdf"), w = 2, h = 3)
p1
p2
p3
p4
dev.off()

pdf(paste0(wd, "/QC-ATAC-Statistics_FRIP.pdf"), w = 2.5, h = 2.5)
p5
p6
dev.off()

#plotPDF(p1, p2, p3, p4, p5, p6, name = "QC-ATAC-Statistics.pdf", ArchRProj = atac_proj, addDOC = FALSE, width = 4, h = 4)

#################################################
# QC plots for combined multiome metrics
#################################################

plot.name <- "QC-Multiome-Statistics"
dir.name <- paste0(wd, "/", plot.name)
dir.create(dir.name, showWarnings = FALSE, recursive = TRUE)

metadata <- atac_proj@cellColData %>% as.data.frame()
readr::write_tsv(metadata, file = paste0(dir.name, "/metadata.tsv"))

samples <- unique(metadata$Sample)

pdf(paste0(plotDir, "/", plot.name, ".pdf"), w = 5, h = 5)
for (sample in samples) {
  df <- filter(metadata, Sample == sample)
  p <- ggplot(df, aes(x = nFrags, y = Gex_nUMI)) +
    ggpointdensity::geom_pointdensity() +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                       breaks=c(1000, 10000, 100000), 
                       labels = expression(10^3, 10^4, 10^5)
    ) +
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                       breaks=c(1000, 10000, 100000), 
                       labels = expression(10^3, 10^4, 10^5)
    ) +
    theme_BOR(border=TRUE) +
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"),
          aspect.ratio = 1,
          legend.position = "none",
          plot.title = element_text(hjust = 0)) +
    xlab("Number of ATAC fragments") +
    ylab("Number of UMIs") +
    ggtitle(sample, 
            subtitle = paste0(sprintf("Median ATAC Fragments = %s\n", median(df$nFrags)),
                              sprintf("Median RNA UMIs = %s", median(df$Gex_nUMI)))) 
  print(p)
}
dev.off()

##########################################################################################
# Dot plots of RNA 
##########################################################################################

# Dot plot of marker genes 

# Markers for identifying fine class of cells
featureSets <- list(
  "Epithelial" = c("EPCAM", "CDH1"), # Epithelial
  "Airway" = c("SOX2", "SCGB3A2"), # Airway
  "PNEC" = c("P3H2", "NRXN1", "ASCL1", "NEUROD1", "GHRL", "GRP"), #PNEC
  "Ciliated" = c("FOXJ1", "DNAH12"), #Ciliated
  "Alveolar" = c("SFTPB", "SOX9", "ETV5", "AGER", "SFTPC"), # Alveolar
  "Cycling" = c("TOP2A", "MKI67"),
  
  "Endothelial" = c("PECAM1", "CLDN5"),
  "gCaps" = c("CA4"),
  "Aerocytes" = c("EDNRB", "S100A3"),
  "Arteries" = c("GJA4", "DKK2", "FBLN5"),
  "Venous" = c("ACKR1", "PLVAP"),
  "Lymphatic" = c("CCL21", "PROX1"),
  
  "Structural" = c("COL1A1", "COL1A2"),
  "Airway" = c("PDGFRA"),
  "AlvF" = c("MEOX2", "WNT2"),
  "AdvF" = c( "PAMR1", "SERPINF1"),
  "MyoF" = c("DACH2", "EYA4"),
  "Contractile" = c("ACTA2", "TAGLN"),
  "AiwaySM" = c("HHIP", "HPSE2"),
  "Vascular" = c("PDGFRB"),
  "Peri" = c("TRPC6", "SULT1E1", "LRRTM4", "RBFOX1", "PAG1"),
  "VascSM" = c("SLIT3", "ELN",  "CSMD1", "PLN", "ITGA11"),
  "Chondrocytes" = c("ACAN", "COL2A1"),
  "Schwann" = c("CDH19", "MPZ"),
  "Mesothelial" = c("C3", "MSLN"),
  
  "Immune" = c("PTPRC"),
  "Mono1" = c("AQP9", "VCAN", "LYZ"),
  "Dc" = c("HLA-DRA", "HLA-DQA1"),
  "IM" = c("MRC1", "MSR1", "C1QB"),
  "Neut" = c("MPO", "LTF"),
  "NK" = c("NKG7", "CCL4"),
  "T_cells" = c("IL7R", "CAMK4"),
  "B_cells" = c("EBF1", "IGHM", "MS4A1"),
  "Plasma" = c("IRF4", "JCHAIN"),
  "MPP" = c("CD34", "GATA2")
) %>% unlist() %>% unname() %>% unique()

FclustOrder <- fineOrder %>% unname() %>% unlist()

Idents(rna_proj) <- "FineNamedClust"
Idents(rna_proj) <- factor(Idents(rna_proj), levels = FclustOrder)

#DefaultAssay(rna_proj) <- "SCT"

pdf(paste0(plotDir, "/RNA_FineNamedClust_markers_dot_plot_lda.pdf"), width = 10, height = 12)
DotPlot(rna_proj, 
        features = rev(featureSets)) + 
  scale_colour_gradientn(colours = cmaps_BOR$sunrise) + #BuenColors::jdb_palette("flame_flame")
  coord_flip() + 
  theme(
    #axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    #legend.title = element_text(size = 5),
    #legend.text = element_text(size = 5),
    axis.title = element_blank()
    ) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dev.off()

# Epithelial only markers
compartment <- "Epithelial"
subFeatureSets <- list(
  "Epithelial" = c("EPCAM", "CDH1"), # Epithelial
  "Airway" = c("SOX2", "SCGB3A2"), # Airway
  "PNEC" = c("P3H2", "NRXN1", "ASCL1", "NEUROD1", "RFX6" ,"GHRL", "GRP"), #PNEC
  "Ciliated" = c("FOXJ1", "DNAH12"), #Ciliated
  "Alveolar" = c("SFTPB", "SFTPC", "SFTPD", "ABCA3", "LAMP3", "SOX9", "ETV5", "AGER"), # Alveolar
  "Cycling" = c("TOP2A", "MKI67")
) %>% unname() %>% unlist()

#"Mesothelial" = c("C3", "MSLN")

subClusters <- c("APr1","APr2","APr3","APr4","AT1l","AT2l","eCili","EpiC","lCili","PNEC1","PNEC2","PNEC3","TiP1","TiP2")

pdf(paste0(plotDir, "/RNA_FineNamedClust_markers_dot_plot_", compartment ,".pdf"), 
    width = length(subFeatureSets) * 0.25, 
    height = length(subClusters) * 0.25)
DotPlot(rna_proj, 
        features = subFeatureSets,
        idents = subClusters) + 
  scale_colour_gradientn(colours = cmaps_BOR$sunrise) + #BuenColors::jdb_palette("flame_flame")
  theme(
    #axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.title = element_text(size = 5),
    #legend.text = element_text(size = 5),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dev.off()

# Endothelial only markers
compartment <- "Endothelial"
subFeatureSets <- list(
  "Endothelial" = c("PECAM1", "CLDN5"),
  "gCaps" = c("CA4", "IL7R"),
  "Aerocytes" = c("GRID1","EDNRB", "FOXF1"),
  "Arteries" = c("GJA4", "DKK2", "FBLN5"),
  "Venous" = c("ACKR1", "HDAC9", "PLVAP"),
  "Lymphatic" = c("CCL21", "PROX1")
) %>% unname() %>% unlist()

#rna_proj[,which(rna_proj$compartment == compartment)]$FineNamedClust %>% table()

subClusters <- c("eAero", "eArtr", "egCap", "lAero", "lArtr", "lgCap",  "Lymp", "Veno")

pdf(paste0(plotDir, "/RNA_FineNamedClust_markers_dot_plot_", compartment ,".pdf"), 
    width = length(subFeatureSets) * 0.3, 
    height = length(subClusters) * 0.3)
DotPlot(rna_proj, 
        features = subFeatureSets,
        idents = subClusters) + 
  scale_colour_gradientn(colours = cmaps_BOR$sunrise) + #BuenColors::jdb_palette("flame_flame")
  theme(
    #axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.title = element_text(size = 5),
    #legend.text = element_text(size = 5),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dev.off()

# Stromal only markers
compartment <- "Stromal"
subFeatureSets <- list(
  "Structural" = c("COL1A1", "COL1A2"),
  "Airway" = c("PDGFRA"),
  "AlvF" = c("MEOX2", "WNT2"),
  "AdvF" = c( "PAMR1", "SERPINF1"),
  "MyoF" = c("DACH2", "EYA4"),
  "Contractile" = c("ACTA2", "TAGLN"),
  "AiwaySM" = c("HHIP", "HPSE2"),
  "Vascular" = c("PDGFRB"),
  "Peri" = c("TRPC6", "SULT1E1", "LRRTM4", "PAG1"),
  "VascSM" = c("SLIT3", "ELN",  "CSMD1", "PLN", "ITGA11"),
  "Chondrocytes" = c("ACAN", "COL2A1"),
  "Schwann" = c("CDH19", "MPZ"),
  "Mesothelial" = c("C3", "MSLN")
) %>% unname() %>% unlist()

subClusters <- rna_proj@meta.data %>% 
  dplyr::filter(compartment == "Stromal") %>% 
  dplyr::select(FineNamedClust) %>% unique() %>% unname() %>% unlist()

subClusters <- c(subClusters, "Meso", "Schw")

pdf(paste0(plotDir, "/RNA_FineNamedClust_markers_dot_plot_", compartment ,".pdf"), 
    width = length(subFeatureSets) * 0.25, 
    height = length(subClusters) * 0.25)
DotPlot(rna_proj, 
        features = subFeatureSets,
        idents = subClusters) + 
  scale_colour_gradientn(colours = cmaps_BOR$sunrise) + #BuenColors::jdb_palette("flame_flame")
  theme(
    #axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.title = element_text(size = 5),
    #legend.text = element_text(size = 5),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dev.off()

# Immune only markers
compartment <- "Immune"
subFeatureSets <- list(
  "Immune" = c("PTPRC"),
  "Mono1" = c("AQP9", "VCAN", "LYZ"),
  "Dc" = c("HLA-DRA", "HLA-DQA1"),
  "IM" = c("MRC1", "MSR1", "C1QB"),
  "Neut" = c("MPO", "LTF"),
  "NK" = c("NKG7", "CCL4"),
  "T_cells" = c("IL7R", "CAMK4"),
  "B_cells" = c("EBF1", "IGHM", "MS4A1"),
  "Plasma" = c("IRF4", "JCHAIN"),
  "MPP" = c("CD34", "GATA2")
) %>% unname() %>% unlist()

subClusters <- rna_proj@meta.data %>% 
  dplyr::filter(compartment == compartment) %>% 
  dplyr::select(FineNamedClust) %>% unique() %>% unname() %>% unlist()

pdf(paste0(plotDir, "/RNA_FineNamedClust_markers_dot_plot_", compartment ,".pdf"), 
    width = length(subFeatureSets) * 0.25, 
    height = length(subClusters) * 0.3)
DotPlot(rna_proj, 
        features = subFeatureSets,
        idents = subClusters) + 
  scale_colour_gradientn(colours = cmaps_BOR$sunrise) + #BuenColors::jdb_palette("flame_flame")
  theme(
    #axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.title = element_text(size = 5),
    #legend.text = element_text(size = 5),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dev.off()

##########################################################################################
# Identifying Marker Peaks FineNamedClust
##########################################################################################

# BroadNamedClust First:
# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
  ArchRProj = atac_proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "FineNamedClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks[,fineOrder], 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  nLabel = 1, # It still seems like there's not actually a way to NOT plot any labels
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  transpose = FALSE
)

pdf(paste0(wd, "/MarkerPeak-Heatmap-CustomNamedClust.pdf"), width = 10, height = 15)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
dev.off()

saveRDS(markerPeaks, paste0(wd, "/MarkerPeaks.rds"))

##########################################################################################
# Motif Enrichments
##########################################################################################

#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = atac_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

saveRDS(enrichMotifs, paste0(wd, "/enrichMotifs.rds"))
#enrichMotifs <- readRDS(paste0(wd, "/enrichMotifs.rds"))

# Rename motifs for more aesthetic plotting:
#rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichMotifs[,FclustOrder], 
  n=5, 
  #clusterCols = FALSE, # This is currently bugged and does not work
  transpose=FALSE, 
  cutOff=log10pCut
)

#draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
#plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-FineNamedClust", width=8, height=12, ArchRProj=atac_proj, addDOC=FALSE)

# Heatmap of motif enrichments
plot_mat <- plotEnrichHeatmap(enrichMotifs[,fineOrder], n=5, transpose=FALSE, 
                              cutOff=log10pCut, returnMatrix=TRUE)

plot_mat <- prettyOrderMat(plot_mat[,fineOrder], clusterCols=FALSE, cutOff=1)$mat

pdf(paste0(plotDir, "/MarkerPeak-MotifEnriched-Heatmap-FineNamedClust.pdf"), width=13, height=18)
fontsize <- 8
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=fineOrder,col=list(atac_cluster=FineNamedClustCmap), 
                        show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  top_annotation = ta,
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

# List of TF motifs used for analysis
tfs <- enrichMotifs@NAMES
tfs[tfs %ni% rownames(rna_proj)] <- c("AC0235093", # Fix names 
                                      "AC1386961",
                                      "MEF2B",
                                      "SOHLH2",
                                      "ZNF875",
                                      "NKX1-1", "NKX1-2", "NKX2-1", "NKX2-2", 
                                      "NKX2-3", "NKX2-4", "NKX2-5", "NKX2-6", 
                                      "NKX2-8", "NKX3-1", "NKX3-2", "NKX6-1", 
                                      "NKX6-2", "NKX6-3",
                                      "T", "ZSCAN5", "ZZZ3"
                                      )
tfs <- unique(tfs)

# Filter for TFs that are expressed in that cell type
avg_expr <- AverageExpression(rna_proj, 
                              assays = "RNA", 
                              group.by = "FineNamedClust", 
                              features = tfs) %>% as.data.frame()

tfs <- rownames(avg_expr)
expressed_tfs <- tfs[rowSums(avg_expr) > 1]

log10pCut <- 20

# Heatmap of motif enrichments
plot_mat <- plotEnrichHeatmap(enrichMotifs[,fineOrder], n=5, transpose=FALSE, 
                              cutOff=log10pCut, returnMatrix=TRUE)

plot_mat <- prettyOrderMat(plot_mat[,fineOrder], clusterCols=FALSE, cutOff=1)$mat

# Clean up the TF names
rownames(plot_mat) <- sub(" .*", "", rownames(plot_mat))

rownames(plot_mat)[rownames(plot_mat) %ni% expressed_tfs] <- c("HNF1A", "CTCFL", "NKX2-6", 
                                                               "NKX2-5", "FOXA3", "FOXB2", 
                                                               "FOXB1", "FOXE3", "FOXI2", 
                                                               "NKX2-1", "GATA1", "TAL2", 
                                                               "ETV2", "SOX30", "TWIST1", 
                                                               "MSC", "ATOH1", "EBF3", "SPIC")

plot_mat_df <- plot_mat %>% as.data.frame()
plot_mat_df <- dplyr::filter(plot_mat_df, rownames(plot_mat_df) %in% expressed_tfs)

pdf(paste0(plotDir, "/MarkerPeak-Expressed-MotifEnriched-Heatmap-FineNamedClust.pdf"), width=13, height=18)
fontsize <- 8
ht_opt$simple_anno_size <- unit(0.25, "cm")
ta <- HeatmapAnnotation(atac_cluster=fineOrder,col=list(atac_cluster=FineNamedClustCmap), 
                        show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
hm <- BORHeatmap(
  plot_mat_df, 
  limits=c(0.0,100.0), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  top_annotation = ta,
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

#Motif enrichment in marker peaks per compartment####
compartments <- atac_proj$compartment %>% unique()

compartments <- c("Stromal")

for (compartment in compartments) {
  subProj <- loadArchRProject(sprintf("/oak/stanford/groups/wjg/skim/projects/LDA/1b_atac_preprocess/subclustering/%s", compartment))

  # Identify Marker Peaks while controling for TSS and Depth Biases
  markerPeaks <- getMarkerFeatures(
    ArchRProj = subProj, 
    useMatrix = "PeakMatrix", 
    groupBy = "FineNamedClust",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  #Visualize Markers as a heatmap
  heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markerPeaks[,fineOrder[[compartment]]], 
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    nLabel = 1, # It still seems like there's not actually a way to NOT plot any labels
    binaryClusterRows = TRUE,
    clusterCols = FALSE,
    transpose = FALSE
  )
  pdf(paste0(plotDir, "/MarkerPeak-Heatmap-", compartment, ".pdf"), width = 5, height = 8)
  draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
  dev.off()
  
  saveRDS(markerPeaks, paste0(plotDir, "/MarkerPeak-", compartment, ".rds"))
  
  #Identify Motif Enrichments
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = subProj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  
  saveRDS(enrichMotifs, paste0(plotDir, "/enrichMotifs-", compartment, ".rds"))
  
  # Subset to clusters that have at least some enrichment
  log10pCut <- 10
  
  #ArchR Heatmap
  heatmapEM <- plotEnrichHeatmap(
    enrichMotifs[,fineOrder[[compartment]]], 
    n=5, 
    #clusterCols = FALSE, # This is currently bugged and does not work
    transpose=FALSE, 
    cutOff=log10pCut
  )
  
  # Heatmap of motif enrichments
  plot_mat <- plotEnrichHeatmap(enrichMotifs[,fineOrder[[compartment]]], n=3, transpose=FALSE, 
                                cutOff=log10pCut, returnMatrix=TRUE)
  
  plot_mat <- prettyOrderMat(plot_mat[,fineOrder[[compartment]]], clusterCols=FALSE, cutOff=1)$mat
  
  pdf(paste0(plotDir, "/MarkerPeak-MotifEnriched-Heatmap-", compartment,".pdf"), width=10, height=14)
  fontsize <- 8
  ht_opt$simple_anno_size <- unit(0.25, "cm")
  ta <- HeatmapAnnotation(atac_cluster=fineOrder[[compartment]],col=list(atac_cluster=FineNamedClustCmap), 
                          show_legend=c(atac_cluster=FALSE), show_annotation_name = c(atac_cluster=FALSE))
  hm <- BORHeatmap(
    plot_mat, 
    limits=c(0.0,100.0), 
    clusterCols=FALSE, clusterRows=FALSE,
    labelCols=TRUE, labelRows=TRUE,
    dataColors = cmaps_BOR$comet,
    top_annotation = ta,
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
}


#################################################
# Dot plots of cluster colors
#################################################

# RNA cluster colors
df <- data.frame(clusters=factor(LNrnaOrder, ordered=TRUE, levels=LNrnaOrder), 
    y=rep(1, length(NrnaOrder)))

pdf(paste0(plotDir, "/scRNA_cluster_colors.pdf"), width=10, height=2)
p <- (ggplot(df, aes(x=clusters, y=y, color=clusters))
  + geom_point(size=10)
  + scale_color_manual(values=rnaLabelClustCmap)
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
    panel.grid.minor= element_blank(), 
    plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
    legend.position = "none", # Remove legend
    axis.text.x = element_text(angle = 90, hjust = 1)) 
)
p
dev.off()


# ATAC cluster colors
df <- data.frame(clusters=factor(LNatacOrder, ordered=TRUE, levels=LNatacOrder), 
    y=rep(1, length(NrnaOrder)))

pdf(paste0(plotDir, "/scATAC_cluster_colors.pdf"), width=10, height=2)
p <- (ggplot(df, aes(x=clusters, y=y, color=clusters))
  + geom_point(size=10)
  + scale_color_manual(values=atacLabelClustCmap)
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
    panel.grid.minor= element_blank(), 
    plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
    legend.position = "none", # Remove legend
    axis.text.x = element_text(angle = 90, hjust = 1)) 
)
p
dev.off()


#################################################
# Violin plots of QC metrics
#################################################

rna_ccd <- rna_proj@meta.data

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

# scRNA nUMIs / cell
p <- (
  ggplot(rna_ccd, aes(x=Sample, y=nCount_RNA, fill=Sample))
  + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
  + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
  + scale_color_manual(values=sample_cmap, limits=names(sample_cmap), name="Sample", na.value="grey")
  + scale_fill_manual(values=sample_cmap)
  + guides(fill=guide_legend(title="Sample"), 
           colour=guide_legend(override.aes = list(size=5)))
  + xlab("")
  + ylab("Number of UMIs per cell")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
  + scale_y_continuous(limits=c(0,25000), expand = c(0, 0))
)

pdf(paste0(plotDir, "/Violin_scRNA_nUMIs.pdf"), width=10, height=4)
p
dev.off()


# scRNA nGenes / cell
p <- (
  ggplot(rna_ccd, aes(x=Sample, y=nFeature_RNA, fill=Sample))
  + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
  + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
  + scale_color_manual(values=sample_cmap, limits=names(sample_cmap), name="Sample", na.value="grey")
  + scale_fill_manual(values=sample_cmap)
  + guides(fill=guide_legend(title="Sample"), 
           colour=guide_legend(override.aes = list(size=5)))
  + xlab("")
  + ylab("Number of genes per cell")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
  + scale_y_continuous(limits=c(0,8000), expand = c(0, 0))
)

pdf(paste0(plotDir, "/Violin_scRNA_nGenes_per_cell.pdf"), width=10, height=4)
p
dev.off()


# scRNA pctMito / cell
p <- (
  ggplot(rna_ccd, aes(x=Sample, y=percent.mt, fill=Sample))
  + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
  + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
  + scale_color_manual(values=sample_cmap, limits=names(sample_cmap), name="Sample", na.value="grey")
  + scale_fill_manual(values=sample_cmap)
  + guides(fill=guide_legend(title="Sample"), 
           colour=guide_legend(override.aes = list(size=5)))
  + xlab("")
  + ylab("Percent Mitochondrial Reads")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
  + scale_y_continuous(limits=c(0,50), expand = c(0, 0))
)

pdf(paste0(plotDir, "/Violin_scRNA_pctMito_per_cell.pdf"), width=10, height=4)
p
dev.off()

atac_ccd <- atac_proj@cellColData %>% as.data.frame()

# scATAC TSS / cell
p <- (
  ggplot(atac_ccd, aes(x=Sample, y=TSSEnrichment, fill=Sample))
  + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
  + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
  + scale_color_manual(values=sample_cmap, limits=names(sample_cmap), name="Sample", na.value="grey")
  + scale_fill_manual(values=sample_cmap)
  + guides(fill=guide_legend(title="Sample"), 
           colour=guide_legend(override.aes = list(size=5)))
  + xlab("")
  + ylab("TSS Enrichment")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
  + scale_y_continuous(limits=c(0,25), expand = c(0, 0))
)

pdf(paste0(plotDir, "/Violin_scATAC_TSS_per_cell.pdf"), width=10, height=4)
p
dev.off()

# scATAC log10 nFrags / cell
p <- (
  ggplot(atac_ccd, aes(x=Sample, y=log10nFrags, fill=Sample))
  + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
  + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
  + scale_color_manual(values=sample_cmap, limits=names(sample_cmap), name="Sample", na.value="grey")
  + scale_fill_manual(values=sample_cmap)
  + guides(fill=guide_legend(title="Sample"), 
           colour=guide_legend(override.aes = list(size=5)))
  + xlab("")
  + ylab("log10 nFrags")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
  + scale_y_continuous(limits=c(0,5), expand = c(0, 0))
)

pdf(paste0(plotDir, "/Violin_scATAC_log10nFrags_per_cell.pdf"), width=10, height=4)
p
dev.off()

#################################################
# Compare cell type accuracy with integration
#################################################

saveArchRProject(atac_proj, outputDirectory = paste0(wd, "/integrated_lda"))

# Unconstrained integration between RNA FineNamedClust and ATAC FineNamedClust
# Compares annotation consistency in terms of gene scores and gene expression
atac_proj <- loadArchRProject(paste0(wd, "/integrated_lda"))
atac_proj <- addGeneIntegrationMatrix(
  ArchRProj = atac_proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "LSI_ATAC",
  seRNA = rna_proj,
  addToArrow = FALSE,
  groupATAC = "FineNamedClust",
  groupRNA = "FineNamedClust",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(atac_proj$FineNamedClust, atac_proj$predictedGroup_Un))
newOrder <- unlist(unname(fineOrder))
cM <- cM[newOrder[which(newOrder %in% rownames(cM))], newOrder[which(newOrder %in% colnames(cM))]]

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}
cM <- new_cM

pdf(paste0(plotDir, "/ATAC-RNA-integration-cM-heatmap.pdf"), width=10, height=10)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  cM, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  row_names_side = "left",
  width = ncol(cM)*unit(0.5, "cm"),
  height = nrow(cM)*unit(0.5, "cm"),
  border_gp=gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()

# Compare ATAC clusters to RNA clustering by integration
atac_proj <- addGeneIntegrationMatrix(
  ArchRProj = atac_proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "LSI_ATAC",
  seRNA = rna_proj,
  addToArrow = FALSE,
  groupATAC = "Clusters_ATAC",
  groupRNA = "FineNamedClust",
  nameCell = "RNA_paired_cell",
  nameGroup = "FineClust_RNA",
  nameScore = "predictedScore"
)

cM <- as.matrix(confusionMatrix(atac_proj$Clusters_ATAC, atac_proj$FineClust_RNA))
#newOrder <- unlist(unname(fineOrder))
#cM <- cM[newOrder[which(newOrder %in% rownames(cM))], newOrder[which(newOrder %in% colnames(cM))]]

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}
cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

pdf(paste0(plotDir, "/ATAC-RNA-integration-FineNamedClust-To-ATAC_Clusters-cM-heatmap.pdf"), width=10, height=10)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  cM, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  row_names_side = "left",
  width = ncol(cM)*unit(0.5, "cm"),
  height = nrow(cM)*unit(0.5, "cm"),
  border_gp=gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()





