#!/bin/bash

# This script calls LongReadSum to generate a quality check analysis summary report of the BAM file following the read filter of the artifact-based filtering step 

# Set the following parameters in the command: 
# 1: Sample Name 


# Example command:   ./09_map_qc_LRS.sh sample_name active
# (active, inactive, or ORF2)


sample_name=$1
L1_ref_type=$2

cd ../$sample_name/d_LINE_quantification/$L1_ref_type/read_filter/

# Remember to change the paths to your specific path!!

BAM_INPUT=../${sample_name}/d_LINE_quantification/${L1_ref_type}/read_filter/${sample_name}_read_filter_passed.sorted_position.bam

python /mnt/isilon/wang_lab/karly/softwares/LongReadSum bam -i $BAM_INPUT -o $sample_name


