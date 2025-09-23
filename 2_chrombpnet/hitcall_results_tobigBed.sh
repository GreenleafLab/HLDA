#!/bin/bash

############################################
# Hit calling results convert to bigBed for visualization
############################################

# Cluster names for hit calling
clusters=(APr1 APr2 APr3 APr4 PNEC1 PNEC2 PNEC3 eCili lCili TiP1 TiP2 AT2l AT1l EpiC Meso egCap lgCap eAero lAero eArtr lArtr Veno Lymp eAlvF lAlvF AdvF MyoF aSMC ePeri lPeri vSMC1 vSMC2 Chdr Mono1 Mono2 IM Dc MPP Tc NK Bc Plasma Neut) #removed Schw since no markerpeaks

# directories
modisco_dir=${wd}/hit_calling/all

output_dir=${wd}/hit_calling/all
# Make directory if not already present
if [ ! -d ${output_dir} ]; then mkdir ${output_dir}; fi

log_dir=${output_dir}/log_dir
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}" # for each cluster in predefined list
do
    jobname=${cluster}_toBB
    out=${log_dir}/${cluster}.toBB.out
    err=${log_dir}/${cluster}.toBB.err
    out_prefix=${output_dir}/${cluster}

    # if the model directory does not exist already
    if [ ! -d ${outdir} ]
    then # extra gpu options with high capacity (GPU_SKU:A40)|(GPU_SKU:V100_SXM2)|
        echo "Hit calling for ${cluster}..."
        sbatch --partition=nih_s10 --gres=gpu:1 --account=wjg --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 64GB --time 1:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "ml cuda/11.2.2_460.32.03 cuda/11.6.0_510.39.01_cudNN_8.7.0.84; finemo call-hits -M pp -r ${output_dir}/${cluster}.npz -m ${modisco_dir}/${cluster}_modisco.h5 -p ${peaks_dir}/${cluster}.bed -o ${outdir}"
        sleep 0.2s
    else
        echo "Hit calling for ${cluster} already exists. Skipping..."
    fi

done

