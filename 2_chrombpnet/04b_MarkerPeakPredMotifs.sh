#!/bin/bash

############################################
# Predict contribution bigwigs for marker peaks
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Cluster names, splits for model
clusters=(APr1 APr2 APr3 APr4 PNEC1 PNEC2 PNEC3 eCili lCili TiP1 TiP2 AT2l AT1l EpiC Meso egCap lgCap eAero lAero eArtr lArtr Veno Lymp eAlvF lAlvF AdvF MyoF aSMC ePeri lPeri vSMC1 vSMC2 Chdr Mono1 Mono2 IM Dc Neut MPP Tc NK Bc Plasma) # removed Schw since no markerpeaks
splits=(fold_2 fold_3 fold_4) # fold_2 fold_3 fold_4

pred_dir=${wd}/markerpeak_motifs_pred/contribution_scores
# Make directory if not already present
if [ ! -d ${pred_dir} ]; then mkdir ${pred_dir}; fi

log_dir=${wd}/markerpeak_motifs_pred/log_dir
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi


for cluster in "${clusters[@]}" # for each cluster in predefined list
do
    for split in "${splits[@]}" # for each split
    do
        regions=${wd}/markerpeak_motifs_pred/peaks/${cluster}.bed
        modelfile=${models_dir}/${cluster}.${split}/models/chrombpnet_nobias.h5 # Grab full model path for peaks
        
        output_dir=${pred_dir}/${cluster}.${split}
        
        jobname=pred_bw_${cluster}_${split}
        output_prefix=${output_dir}/${cluster}.${split}
        out=${log_dir}/${cluster}.${split}.out
        err=${log_dir}/${cluster}.${split}.err

        # if the model directory does not exist already
        if [ ! -d ${output_dir} ]
        then
            mkdir ${output_dir}
            echo "Predicting bigwigs for ${cluster} with ${split}..."
            sbatch -p gpu --gpus 1 --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 100GB --time 48:00:00 --no-requeue \
            --job-name=${jobname} --output=${out} --error=${err} \
            --wrap "ml cudnn/8.1; ml cuda/11.2.0; chrombpnet contribs_bw -m ${modelfile} -r ${regions} -g ${genome} -c ${chromsize} -op ${output_prefix}"
            sleep 1s
        else
            echo "Contribution scores for ${cluster} with ${split} already exists. Skipping..."
        fi
    done
done
