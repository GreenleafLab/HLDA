#!/bin/bash

############################################
# De novo motif discovery with TF Modisco
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate new_chrombpnet

# Cluster names, splits for model
clusters=(APr1 APr2 APr3 APr4 AT2l AT1l Bc Chdr Dc eAero eCili egCap EpiC IM lCili lgCap Lymp Meso Mono1 Mono2 MPP Neut Plasma PNEC1 PNEC2 PNEC3 Schw Tc TiP1 Veno vSMC1 vSMC2 TiP2 lAero eArtr lArtr eAlvF lAlvF AdvF MyoF ePeri aSMC lPeri NK)

# Directory with averaged contribution scores for each cluster
base_dir=${wd}/fullpeak_pred/contribution_scores

# Output directory
modisco_dir="${wd}/fullpeak_pred/modisco_outputs"
# Make directory for if not already
if [ ! -d ${modisco_dir} ]; then mkdir ${modisco_dir}; fi

# Log directory
log_dir="${modisco_dir}/log_dir_modisco"
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi


for cluster in "${clusters[@]}"
do            
    h5py=${base_dir}/${cluster}_average_counts_scores.h5
    jobname=${cluster}_modisco
    out=${log_dir}/${cluster}.out
    err=${log_dir}/${cluster}.err

    output_prefix=${modisco_dir}/${cluster}
    outfile=${modisco_dir}/${cluster}_modisco.h5

    # if the model directory does not exist already
    if [ ! -f ${outfile} ]
    then
        echo "Motif discovery with TF Modisco for ${cluster}..."
        sbatch -p wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 100GB \
        --time 48:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "chrombpnet modisco_motifs -i ${h5py} -n 1000000 -op ${output_prefix} -v"
        sleep 10s
    else
        echo "TF Modisco results for ${cluster} already exists. Skipping..."
    fi
done
