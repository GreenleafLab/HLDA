#!/bin/bash

############################################
# Create BigWigs from averaged contribution scores
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate new_chrombpnet

# Cluster names, splits for model
clusters=(egCap lgCap eAero lAero Lymp eArtr lArtr Veno PNEC1 PNEC2 PNEC3) # 

# Python script for hdf5 to bigwig
python_path="/oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/importance_hdf5_to_bigwig.py"

# Directory with averaged contribution scores for each cluster
base_dir=${wd}/pos_ctrl_pred/contribution_scores

# Log directory
log_dir=${base_dir}/log_dir_BigWigAverage
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}"
do            
    filename=${cluster}
    h5py=${base_dir}/${filename}_average_counts_scores.h5
    bedfile="/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/pos_ctrl_pred/PROX1_enh_peak.bed"

    jobname=${filename}_bigwig
    out=${log_dir}/${filename}.out
    err=${log_dir}/${filename}.err

    out_prefix=${base_dir}/${filename}_average_counts_scores
    outfile=${base_dir}/${filename}_average_counts_scores.bw

    # if the bigwig file does not exist already
    if [ ! -f ${outfile} ]
    then
        echo "Contribution score bigwig for ${filename}..."
        sbatch -p wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 1 --mem 8GB \
        --time 1:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "python ${python_path} -h5 ${h5py} -r ${bedfile} -c ${chromsize} -op ${out_prefix}"
        sleep 0.2s    
    else
        echo "Contribution score bigwig for ${filename} already exists. Skipping..."
    fi
done
