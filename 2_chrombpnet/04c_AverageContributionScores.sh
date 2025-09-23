#!/bin/bash

############################################
# Predict contribution bigwigs for marker peaks
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Dependency
# Requires average_h5.py to perform the per cell type averaging
avg_script_path="/oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/average_h5.py"

# Cluster names
clusters=(APr1 APr2 APr3 APr4 PNEC1 PNEC2 PNEC3 eCili lCili TiP1 TiP2 AT2l AT1l EpiC Meso egCap lgCap eAero lAero eArtr lArtr Veno Lymp eAlvF lAlvF AdvF MyoF aSMC ePeri lPeri vSMC1 vSMC2 Chdr Mono1 Mono2 IM Dc MPP Tc NK Bc Plasma Neut) #removed Schw since no markerpeaks

# Directory where contribution scores are outputed to
base_dir=${wd}/markerpeak_motifs_pred/contribution_scores

log_dir=${wd}/markerpeak_motifs_pred/log_dir_average
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}" # for each cluster in predefined list
do            
    jobname=${cluster}_avg
    out=${log_dir}/${cluster}.out
    err=${log_dir}/${cluster}.err
    output_file=${base_dir}/${cluster}_average_counts_scores.h5

    if [ ! -d ${output_file} ]
    then
        # if the model directory does not exist already
        echo "Averaging contribution scores for ${cluster}..."
        sbatch -p wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 1 --mem 8GB --time 1:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "python ${avg_script_path} ${cluster} ${base_dir}"
        sleep 0.2s
    else
        echo "Averaged contribution score file for ${cluster} already exists. Skipping..."
    fi
done