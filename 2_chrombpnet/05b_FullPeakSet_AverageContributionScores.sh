#!/bin/bash

############################################
# Average contribution scores
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate new_chrombpnet

# Dependency
# Requires average_h5.py to perform the per cell type averaging
avg_script_path="/oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/average_h5.py"

# Cluster names
clusters=(APr1 APr2 APr3 APr4 AT2l AT1l Bc Chdr Dc eAero eCili egCap EpiC IM lCili lgCap Lymp Meso Mono1 Mono2 MPP Neut Plasma PNEC1 PNEC2 PNEC3 Schw Tc TiP1 Veno vSMC1 vSMC2 TiP2 lAero eArtr lArtr eAlvF lAlvF AdvF MyoF ePeri aSMC lPeri NK)

# Directory where contribution scores are outputed to
base_dir=${wd}/fullpeak_pred/contribution_scores

log_dir=${wd}/fullpeak_pred/log_dir_average
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}" # for each cluster in predefined list
do            
    jobname=${cluster}_avg
    out=${log_dir}/${cluster}.out
    err=${log_dir}/${cluster}.err
    output_file=${base_dir}/${cluster}_average_counts_scores.h5

    if [ ! -f ${output_file} ]
    then
        # if the model directory does not exist already
        echo "Averaging contribution scores for ${cluster}..."
        sbatch -p wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 4 --mem 64GB --time 6:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "python ${avg_script_path} ${cluster} ${base_dir}"
        sleep 5s
    else
        echo "Averaged contribution score file for ${cluster} already exists. Skipping..."
    fi
done
