#!/bin/bash

############################################
# Create BigWigs from averaged contribution scores
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate new_chrombpnet

# Cluster names, splits for model
clusters=(APr1 APr2 APr3 APr4 AT2l AT1l Bc Chdr Dc eAero eCili egCap EpiC IM lCili lgCap Lymp Meso Mono1 Mono2 MPP Neut Plasma PNEC1 PNEC2 PNEC3 Schw Tc TiP1 Veno vSMC1 vSMC2 TiP2 lAero eArtr lArtr eAlvF lAlvF AdvF MyoF ePeri aSMC lPeri NK)

# Python script for hdf5 to bigwig
python_path="/oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/importance_hdf5_to_bigwig.py"

# Directory with averaged contribution scores for each cluster
base_dir=${wd}/fullpeak_pred/contribution_scores

# Log directory
log_dir=${base_dir}/log_BigWigAverage
# Make directory if not already present
if [ ! -d ${log_dir} ]; then mkdir ${log_dir}; fi

for cluster in "${clusters[@]}"
do            
    h5py=${base_dir}/${cluster}_average_counts_scores.h5
    bedfile=${wd}/peaks/${cluster}.bed

    jobname=${cluster}_bigwig
    out=${log_dir}/${cluster}.out
    err=${log_dir}/${cluster}.err

    out_prefix=${base_dir}/${cluster}_average_counts_scores
    outfile=${base_dir}/${cluster}_average_counts_scores.bw

    # if the bigwig file does not exist already
    if [ ! -f ${outfile} ]
    then
        echo "Contribution score bigwig for ${cluster}..."
        sbatch -p owners,wjg,sfgf,biochem,crabtree --mail-type=FAIL --mail-user=samkim93@stanford.edu -c 4 --mem 32GB \
        --time 6:00:00 --no-requeue \
        --job-name=${jobname} --output=${out} --error=${err} \
        --wrap "python ${python_path} -h5 ${h5py} -r ${bedfile} -c ${chromsize} -op ${out_prefix}"
        sleep 10s    
    else
        echo "Contribution score bigwig for ${cluster} already exists. Skipping..."
    fi
done
