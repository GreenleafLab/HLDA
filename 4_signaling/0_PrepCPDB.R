#!/usr/bin/env Rscript

##################################################
# Convert seurat object to anndata for cellphonedb
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
wd <- "/oak/stanford/groups/wjg/skim/projects/LDA/4_signaling/data"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

##### Read in data #####
message("Reading in data...")
obj <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/1a_rna_preprocess/lda_v2.rds")
DefaultAssay(obj) <- "RNA"

##### Convert to AnnData #####
sceasy::convertFormat(obj, from="seurat", to="anndata", outFile='lda.h5ad', assay = "RNA")

obj@meta.data$Cell <- rownames(obj@meta.data)
df = obj@meta.data[, c("Cell", "FineNamedClust")]
write.table(df, file = paste0(wd, "/meta.tsv"), sep = "\t", quote = F, row.names = F)

##### Convert to AnnData for vSMC relevant cell types only #####
Idents(obj) <- "FineNamedClust"
subobj <- subset(obj, idents = c("eArtr", "lArtr", "Veno", "vSMC1", "vSMC2"))
subobj$FineNamedClust[which(subobj$FineNamedClust %in% c("eArtr", "lArtr"))] <- "Artr"

sceasy::convertFormat(subobj, from="seurat", to="anndata", outFile='vSMC_vascular.h5ad', assay = "RNA")
subobj@meta.data$Cell <- rownames(subobj@meta.data)
df = subobj@meta.data[, c("Cell", "FineNamedClust")]
write.table(df, file = paste0(wd, "/vSMC_vascular_meta.tsv"), sep = "\t", quote = F, row.names = F)

##### Prepare DE Genes #####
library(DESeq2)

DE_dir <- paste0(wd, "/DE_Genes")
dir.create(DE_dir, showWarnings = FALSE, recursive = TRUE)

# Load prior vSMC DE gene analysis
resOrdered.vSMC <- readRDS("/oak/stanford/groups/wjg/skim/projects/LDA/Figure03_Stromal/DE_testing/resOrdered.rds") %>% as.data.frame()
resOrdered.vSMC$gene <- rownames(resOrdered.vSMC)

# thresholds
l2fc_threshold <- 1.0
padj_threshold <- 0.01

# filter gene lists 
res.df.vSMC <- resOrdered.vSMC %>% dplyr::filter(abs(log2FoldChange) > l2fc_threshold & padj < padj_threshold)

res.df.vSMC$cluster <- "vSMC1"
res.df.vSMC$cluster[which(res.df.vSMC$log2FoldChange > l2fc_threshold & res.df.vSMC$padj < padj_threshold)] <- "vSMC2"


# Perform DE for vascular cell types
# pseudobulk the counts based on sample and cell type
subobj$vascular_name <- subobj$FineNamedClust
subobj$vascular_name[which(subobj$vascular_name %in% c("eArtr", "lArtr"))] <- "Artr"

pseudo_rna <- AggregateExpression(subobj, 
                                  assays = "RNA", 
                                  return.seurat = T, 
                                  group.by = c("Sample", "vascular_name"), slot = "counts")
pseudo_rna$Sample <- gsub("_(.*)", "", Cells(pseudo_rna))
pseudo_rna$FineNamedClust <- gsub(".*_", "", Cells(pseudo_rna)) 
Idents(pseudo_rna) <- "FineNamedClust"

pseudo_rna_sub <- subset(pseudo_rna, ident = c("Artr", "Veno"))
data.use <- GetAssayData(pseudo_rna_sub, "counts", "RNA")

# Add metadata information on sample (batch)
group.info <- data.frame(celltype = as.factor(gsub(".*_", "", unname(colnames(data.use)))),
                         samplename = unname(colnames(data.use)), 
                         sample = gsub("_(.*)", "", unname(colnames(data.use)))
                         )
group.info$celltype <- relevel(group.info$celltype, ref = "Veno")

dds1 <- DESeq2::DESeqDataSetFromMatrix(
  countData = data.use,
  colData = group.info,
  design = ~ celltype + sample
)

# Run DESeq2
dds1 <- DESeq(dds1)
resultsNames(dds1)
resLFC <- lfcShrink(dds1, coef="celltype_Artr_vs_Veno", type="apeglm")
resOrdered <- resLFC[order(resLFC$pvalue),]

saveRDS(dds1, paste0(DE_dir, "/DDS_Artr_vs_Veno.rds"))
saveRDS(resOrdered, paste0(DE_dir, "/resOrdered_Artr_vs_Veno.rds"))
#resOrdered <- readRDS(paste0(DE_dir, "/resOrdered.rds"))

res.df.vasc <- resOrdered %>% as.data.frame() %>% dplyr::filter(abs(log2FoldChange) > l2fc_threshold & padj < padj_threshold)
res.df.vasc$gene <- rownames(res.df.vasc)

res.df.vasc$cluster <- "Veno"
res.df.vasc$cluster[which(res.df.vasc$log2FoldChange > l2fc_threshold & res.df.vasc$padj < padj_threshold)] <- "Artr"

res.df <- rbind(res.df.vasc, res.df.vSMC)
res.df <- res.df[,c("cluster", "gene", "padj", "log2FoldChange", "baseMean")]
write.table(res.df, file =paste0(DE_dir, "/DEGs_vSMC_vascular.tsv"), sep = '\t', quote = F, row.names = F)

##### Prepare microenvironments file #####
# Based on FISH, vSMC1 is colocalized with Artr and vSMC2 is colocalized with Veno

microenv.df <- data.frame(
	cell_type = c("Artr", "vSMC1", "Veno", "vSMC2"),
	microenvironment= c("Artr|vSMC1", "Artr|vSMC1", "Veno|vSMC2", "Veno|vSMC2")
	)

write.table(microenv.df, file =paste0(wd, "/vSMC_microenvironment.tsv"), sep = '\t', quote = F, row.names = F)

##### Prepare active TF file #####
tf.df <- data.frame(
	cell_type = c("Artr", "vSMC1", "Veno", "vSMC2"),
	TF = c(c("a", "a", "a"), c("b","b","b"))
	)





