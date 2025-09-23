#!/bin/bash

############################################
# Create bias model for each fold
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Make bias models directory if it doesn't already exist
if [ ! -d ${bias_dir} ]; then mkdir ${bias_dir}; fi

# Use the cluster after testing the top 5 cell types with the most number of reads
cluster=egCap

for json in ${splits_dir}/*.json # for each split file
do
  split=$(basename ${json} | sed 's/.json//') # grab fold number
	frag=${frags_dir}/FineNamedClust.${cluster}.tsv.gz
	peaks=${peaks_dir}/${cluster}.bed
	nonpeaks=${nonpeaks_dir}/${cluster}.${split}_negatives.bed
  outfile=${bias_dir}/${cluster}.${split}
  mode=ATAC
  parameter=0.5
    
  out=${bias_dir}/${cluster}.${split}.out
  err=${bias_dir}/${cluster}.${split}.err

  # if the bias model has not been made already
  if [ ! -d ${outfile} ]
  then
      mkdir ${outfile}
      echo "Training bias models for ${split} based on ${cluster}..."
      sbatch -p gpu --gpus 1 --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 10 --mem 80GB \
      --constraint="(GPU_SKU:RTX_2080Ti)|(GPU_SKU:A100_PCIE)|(GPU_SKU:A100_SXM4)" \
      --time 1-0 \
      --job-name=biasmodel_${cluster}_${split} --output=${out} --error=${err} \
      --wrap "ml cudnn/8.1; ml cuda/11.2.0; chrombpnet bias pipeline -ifrag ${frag} -d ${mode} -g ${genome} -c ${chromsize} -p ${peaks} -n ${nonpeaks} -fl ${json} -b ${parameter} -o ${outfile} -fp ${split}"
      sleep 0.25s
  else
      echo "Bias model for ${split} based on ${cluster} already exists. Skipping..."
  fi
done