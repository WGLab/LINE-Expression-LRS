#!/bin/bash

# This script calls LongReadSum to generate a quality check analysis summary report of the mapping to the hg38 reference genome 

# Set the following parameters in the command: 
# 1: Sample Name 

# Example command:   ./05_map_qc_LRS.sh sample_name


sample_name=$1

cd ../$sample_name/c_hg38_mapping_LRS/

# Remember to change the paths to your specific path!!
BAM_INPUT=${sample_name}_hg38_mapped.sorted_position.bam

python /mnt/isilon/wang_lab/karly/softwares/LongReadSum bam -i $BAM_INPUT -o $sample_name


