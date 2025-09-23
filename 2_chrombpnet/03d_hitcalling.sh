#!/bin/bash

############################################
# Hit calling
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate finemo

# Cluster names for hit calling
clusters=(egCap lgCap eAero lAero Lymp eArtr lArtr Veno PNEC1 PNEC2 PNEC3) # 

# Peak file used for contribution score calculations
peakfile="/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/pos_ctrl_pred/PROX1_enh_peak.bed"

# directories
modisco_dir=${wd}/fullpeak_pred/modisco_outputs

output_dir="/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/pos_ctrl_pred/hit_calling"
# Make directory if not already present
if [ ! -d ${output_dir} ]; then mkdir ${output_dir}; fi

log_dir=${output_dir}/log_dir
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}" # for each cluster in predefined list
do
    jobname=${cluster}_hitcall
    out=${log_dir}/${cluster}.hitcall.out
    err=${log_dir}/${cluster}.hitcall.err
    outdir=${output_dir}/${cluster}

    # if the model directory does not exist already
    if [ ! -d ${outdir} ]
    then 
        echo "Hit calling for ${cluster}..."
        sbatch -p owners,gpu,wjg --gpus 1 --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 32GB --time 1:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "ml cudnn/8.1; ml cuda/11.2.0; finemo call-hits -M pp -r ${output_dir}/${cluster}.npz -m ${modisco_dir}/${cluster}_modisco.h5 -p ${peakfile} -o ${outdir}"
        sleep 0.2s
    else
        echo "Hit calling for ${cluster} already exists. Skipping..."
    fi

done

