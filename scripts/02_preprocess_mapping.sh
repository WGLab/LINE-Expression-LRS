#!/bin/bash

# This script is part of the data preprocessing and maps the input reads to the Custom LINE Reference Library, removing reads that do not align to a LINE/L1 element 

# Set the following parameters in the command: 
# 1: Sample Name (no_spaces)

# Example command:   ./02_preprocess_mapping.sh sample_name

sample_name=$1       # this is the name of the sample you are analyzing 
echo "$sample_name"


cd ../$sample_name/a_dataset

echo "Mapping to Custom LINE Reference Library ..."

REF_L1_mega=../../references/custom_LINE_reference.fasta

minimap2 -ax map-ont $REF_L1_mega $1_cDNA_1kb.fasta -t 24 > $1_mapped_cDNA_1kb.sam
samtools fasta $1_mapped_cDNA_1kb.sam -F 2308 > $1_mapped_cDNA_1kb.fa

echo "Mapping to Custom LINE Reference Library complete!" 





