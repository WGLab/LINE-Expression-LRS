#!/bin/bash

# This script utilizes RepeatMasker for LINE/L1 detection, post-processing of the output file, and artifact-based filtering for L1 quantification. 

# Set the following parameters in the command: 
# 1: Sample Name 

# Example command:   ./03_L1_detection.sh sample_name



sample_name=$1


output_dir=../$sample_name/b_repeat_masker_process/
input_file=../$sample_name/a_dataset/${sample_name}_mapped_cDNA_1kb.fa

input_fasta=../a_dataset/${sample_name}_cDNA_1kb.fasta


echo "Running RepeatMasker..."
RepeatMasker -pa 6 -dir ${output_dir} -nolow -norna -div 10 -species human -no_is -a -u -xsmall -xm ${input_file}

echo "Finished RepeatMasker!"



# Post-RepeatMasker Filtering by 10% Divergence


cd ../$sample_name/b_repeat_masker_process/
RM_input_file=${sample_name}_mapped_cDNA_1kb.fa.out
output_file="div10_LINEs.out"

readIDs_lines=()

while read -r line; do
    if [[ $line == *"LINE/L1"* ]]; then
        fields=($line)

        if (( $(echo "${fields[1]} <= 10" | bc -l) )); then
            echo "$line" >> "$output_file"
            readIDs_lines+=("${fields[4]}")
        fi
    fi
done < "$RM_input_file"



# Get ReadIDs and Generate the new FASTA file of these ReadIDs
read_id_file="div10_readIDs.txt"

echo "Gathering Final ReadIDs of less than 10% diverged LINE/L1 elements..."

for item in "${readIDs_lines[@]}"; do
    echo "$item" >> "$read_id_file"
done


seqtk subseq $input_fasta div10_readIDs.txt > ${sample_name}_div10.fa

echo "Completed processing RepeatMasker output and generated a new FASTA file of reads with less than 10% diverged LINE/L1 elements!"

