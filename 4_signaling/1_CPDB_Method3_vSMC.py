# CellPhoneDB on human lung development atlas

# Load packages
import pandas as pd
import sys
import os
import numpy as np
import anndata

pd.set_option('display.max_columns', 100)
os.chdir('/oak/stanford/groups/wjg/skim/projects/LDA/final/signaling')
print(sys.version)

# Input files
cpdb_file_path = '/oak/stanford/groups/wjg/skim/resources/cellphonedb/v4.1.0/cellphonedb.zip'
meta_file_path = '/oak/stanford/groups/wjg/skim/projects/LDA/final/signaling/data/vSMC_vascular_meta.tsv'
counts_file_path = '/oak/stanford/groups/wjg/skim/projects/LDA/final/signaling/data/vSMC_vascular.h5ad'
microenvs_file_path = '/oak/stanford/groups/wjg/skim/projects/LDA/final/signaling/microenv/vSMC_microenv.tsv'
degs_file_path = '/oak/stanford/groups/wjg/skim/projects/LDA/final/signaling/data/DEGs_vSMConly.tsv'
out_path = 'results/method3_vSMCs_only'

# Check metadata files
metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata.head(3)

# Check anndata object
adata = anndata.read_h5ad(counts_file_path)
adata.shape

# Check barcode concordance
list(adata.obs.index).sort() == list(metadata["Cell"]).sort()

# DEGs
pd.read_csv(degs_file_path,
            sep = '\t').head(3)

# Run method 3 analysis
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                            # mandatory: CellPhoneDB database zip file.
    meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,                        # mandatory: normalized count matrix.
    degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
    counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
    microenvs_file_path = microenvs_file_path,                  # optional (default: None): defines cells per microenvironment.
    threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
    separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_path = out_path,                                     # Path to save results
    output_suffix = None                                        # Replaces the timestamp in the output files by a user defined string in the  (default: None)
    )

# Plot exploratory plots
import ktplotspy as kpy
import matplotlib.pyplot as plt

# Plot the total number of significant interactions for each pair of cells
kpy.plot_cpdb_heatmap(
        adata = adata,
        pvals = relevant_interactions,
        celltype_key = "FineNamedClust",
        figsize = (5,5),
        title = "Number of significant interactions",
        symmetrical = False,
        degs_analysis = True
    )

plt.savefig("results/method3_vSMCs_only/TotalSignificantInteractions.svg")
