#!/bin/bash

############################################
# Train chrombpnet model for each cluster
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Make directory for chrombpnet models if not already
if [ ! -d ${models_dir} ]; then mkdir ${models_dir}; fi

# Cluster names, splits for model
clusters=(APr1 APr2 APr3 APr4 PNEC1 PNEC2 PNEC3 eCili lCili TiP1 TiP2 AT2l AT1l EpiC Schw Meso egCap lgCap eAero lAero eArtr lArtr Veno Lymp eAlvF lAlvF AdvF MyoF aSMC ePeri lPeri vSMC1 vSMC2 Chdr Mono1 Mono2 IM Dc Neut MPP Tc NK Bc Plasma)
#clusters=(aSMC)
splits=(fold_0 fold_1 fold_2 fold_3 fold_4) #fold_0 fold_1 fold_2 fold_4

for cluster in "${clusters[@]}" # for each cluster in predefined list
do
    peakfile=${peaks_dir}/${cluster}.bed # Grab full bedfile path for peaks

    for split in "${splits[@]}" # for each split
    do
        bias_model_file=${bias_dir}/egCap.fold_0_bias.h5 # Grab full bias model path   ################################Change after bias model training complete
        json=${splits_dir}/${split}.json
        frag=${frags_dir}/FineNamedClust.${cluster}.tsv.gz
        mode=ATAC
        nonpeaks=${nonpeaks_dir}/${cluster}.${split}_negatives.bed

        outfile=${models_dir}/${cluster}.${split}
        out=${outfile}.out
        err=${outfile}.err
        jobname=cbpnet_${cluster}_${split}

        # if the model directory does not exist already
        if [ ! -d ${outfile} ]
        then
            echo "Training chrombpnet of ${cluster} with ${split}..."
            sbatch -p gpu --gpus 1 --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 75GB \
            --constraint="(GPU_SKU:RTX_2080Ti)|(GPU_SKU:A100_PCIE)|(GPU_SKU:A100_SXM4)" --time 2-0 --no-requeue \
            --job-name=${jobname} --output=${out} --error=${err} \
            --wrap "ml cudnn/8.1; ml cuda/11.2.0; chrombpnet pipeline -ifrag ${frag} -d ${mode} -g ${genome} -c ${chromsize} -p ${peakfile} -n ${nonpeaks} -fl ${json} -b ${bias_model_file} -o ${outfile}"
            sleep 1s
        else
            echo "Chrombpnet model for ${split} based on ${cluster} already exists. Skipping..."
        fi
    done
done
