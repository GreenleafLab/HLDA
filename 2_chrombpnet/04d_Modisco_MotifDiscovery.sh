#!/bin/bash

############################################
# De novo motif discovery with TF Modisco
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Cluster names, splits for model
clusters=(APr1) 
#APr2 APr3 APr4 PNEC1 PNEC2 PNEC3 eCili lCili TiP1 TiP2 AT2l AT1l EpiC Meso egCap lgCap eAero lAero eArtr lArtr Veno Lymp eAlvF lAlvF AdvF MyoF aSMC ePeri lPeri vSMC1 vSMC2 Chdr Mono1 Mono2 IM Dc Neut MPP Tc NK Bc Plasma
# removed Schw since no markerpeaks

# Directory with averaged contribution scores for each cluster
base_dir=${wd}/markerpeak_motifs_pred/contribution_scores

# Output directory
modisco_dir=${wd}/markerpeak_motifs_pred/modisco_outputs
# Make directory for if not already
if [ ! -d ${modisco_dir} ]; then mkdir ${modisco_dir}; fi

# Log directory
log_dir=${wd}/markerpeak_motifs_pred/log_dir_modisco
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi


for cluster in "${clusters[@]}"
do            
    h5py=${base_dir}/${cluster}_average_counts_scores.h5
    jobname=${cluster}_modisco
    out=${log_dir}/${cluster}.out
    err=${log_dir}/${cluster}.err

    outfile=${modisco_dir}/${cluster}.modisco_results.h5

    # if the model directory does not exist already
    if [ ! -d ${outfile} ]
    then
        echo "Motif discovery with TF Modisco for ${cluster}..."
        sbatch -p wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 100GB \
        --time 48:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "modisco motifs -i ${h5py} -n 1000000 -o ${outfile} -v"
        sleep 0.2s    
    else
        echo "Directory for TF Modisco for ${cluster} already exists. Skipping..."
    fi
done
