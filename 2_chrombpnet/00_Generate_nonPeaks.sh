#!/bin/bash

############################################
# Create non-peaks for each peak set
############################################

# Source config file
source /oak/stanford/groups/wjg/skim/projects/LDA/scripts/chrombpnet/00_chrombpnet_config.conf

# Activate conda environment
source /oak/stanford/groups/wjg/skim/software/miniconda3/etc/profile.d/conda.sh
conda activate chrombpnet

# Temp
#bedfiles=(/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/peaks/APr2.bed /oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/peaks/PNEC3.bed /oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/peaks/Cili.bed)

# Make nonpeaks directory if it doesn't already exist
if [ ! -d ${nonpeaks_dir} ]; then mkdir ${nonpeaks_dir}; fi

for bf in ${peaks_dir}/*.bed # for each bed file
#for bf in "${bedfiles[@]}" # for each bed file
do
    cluster=$(basename ${bf} | sed 's/.bed//') # grab just the name of each bedfile

    for json in ${splits_dir}/*.json # for each split file
    do
    	split=$(basename ${json} | sed 's/.json//') # grab fold number
    	prefix=${nonpeaks_dir}/${cluster}.${split}

    	echo "Making non-peaks for ${split} of ${cluster}..."
    	out=${nonpeaks_dir}/${cluster}.${split}.out
    	err=${nonpeaks_dir}/${cluster}.${split}.err

    	sbatch -p wjg,sfgf,biochem,crabtree -t 2:00:00 --mem=16G --cpus-per-task=1 --mail-type=FAIL --mail-user=samkim93@stanford.edu\
    	--job-name=nonpeaks_${cluster}_${split} --output=${out} --error=${err} \
    	--wrap "chrombpnet prep nonpeaks \
    	-g ${genome} \
    	-p ${bf} \
    	-c ${chromsize} \
    	-fl ${json} \
    	-br ${blacklist} \
    	-o ${prefix}"
    	#-il 500" # default = 2114 but changed since ArchR peaks are fixed widths
    	sleep 0.25s
    done
done
