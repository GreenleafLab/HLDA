#!/bin/bash

############################################
# Average contribution scores
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
clusters=(egCap lgCap eAero lAero Lymp eArtr lArtr Veno PNEC1 PNEC2 PNEC3) # 

# Directory with contribution scores 
base_dir=${wd}/pos_ctrl_pred/contribution_scores

# For logging outputs and errors
log_dir=${wd}/pos_ctrl_pred/log_dir_average
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}" # for each cluster in predefined list
do
    filename=${cluster}
    jobname=${filename}_avg
    out=${log_dir}/${filename}.out
    err=${log_dir}/${filename}.out
    output_file=${base_dir}/${filename}_average_counts_scores.h5

    if [ ! -d ${output_file} ]
    then
        # if the model directory does not exist already
        echo "Averaging contribution scores for ${filename}..."
        sbatch -p wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 1 --mem 8GB --time 1:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "python ${avg_script_path} ${filename} ${base_dir}"
        sleep 0.2s
    else
        echo "Averaged contribution score file for ${filename} already exists. Skipping..."
    fi
done
