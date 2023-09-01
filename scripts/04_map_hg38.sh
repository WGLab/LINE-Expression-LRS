#!/bin/bash

# This script maps the reads with LINE/L1 elements less than 10% diverged (as identified from RepeatMasker) followed by LongReadSum to generate a summary report of the mapping quality 

# Set the following parameters in the command: 
# 1: Sample Name 

sample_name=$1

# Example command:   ./04_map_hg38.sh sample_name



# Remember to change the paths to your specific path!!

FASTA_INPUT=../$sample_name/b_repeat_masker_process/${sample_name}_div10.fa
REF_SPLICE=../references/gencode.v40.annotation.bed
REF_GENOME38=../references/hg38.fa

cd ../$sample_name/c_hg38_mapping_LRS/

echo "Mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome..."

minimap2 -ax splice --junc-bed $REF_SPLICE -uf --secondary=no -k14 -t 8 $REF_GENOME38 $FASTA_INPUT -o ${sample_name}_hg38_mapped.sam
samtools view -Sb -o ${sample_name}_hg38_mapped.bam ${sample_name}_hg38_mapped.sam
samtools sort ${sample_name}_hg38_mapped.bam -o ${sample_name}_hg38_mapped.sorted_position.bam 
samtools index ${sample_name}_hg38_mapped.sorted_position.bam 

echo "Finished mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome" 

