#!/bin/bash

############################################
# Train chrombpnet model for each cluster
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Cluster names, splits for model
clusters=(egCap lgCap eAero lAero Lymp eArtr lArtr Veno PNEC1 PNEC2 PNEC3) # 
splits=(fold_0 fold_1 fold_2 fold_3 fold_4) #fold_0 fold_1 fold_2 fold_4

pred_dir=${wd}/pos_ctrl_pred/contribution_scores
# Make directory for chrombpnet models if not already
if [ ! -d ${pred_dir} ]; then mkdir ${pred_dir}; fi

log_dir=${wd}/pos_ctrl_pred/log_dir
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

# regions to predict
regions="/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/pos_ctrl_pred/PROX1_enh_peak.bed"

for cluster in "${clusters[@]}" # for each cluster in predefined list
do
    for split in "${splits[@]}" # for each split
    do
        modelfile=${models_dir}/${cluster}.${split}/models/chrombpnet_nobias.h5 # Grab full model path for peaks
        output_dir=${pred_dir}/${cluster}.${split}
        output_prefix=${output_dir}/${cluster}.${split}
        out=${log_dir}/${cluster}.${split}.out
        err=${log_dir}/${cluster}.${split}.err
        jobname=${cluster}_${split}_pred_bw

        # if the model directory does not exist already
        if [ ! -d ${output_dir} ]
        then
            mkdir ${output_dir}
            echo "Predicting contribution scores for ${cluster} with ${split}..."
            sbatch -p owners,gpu --gpus 1 --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 20GB \
            --time 1:00:00 --no-requeue \
            --job-name=${jobname} --output=${out} --error=${err} \
            --wrap "ml cudnn/8.1; ml cuda/11.2.0; chrombpnet contribs_bw -m ${modelfile} -r ${regions} -g ${genome} -c ${chromsize} -op ${output_prefix}"
            sleep 1s
        else
            echo "Contribution scores for ${cluster} with ${split} already exists. Skipping..."
        fi
    done
done
