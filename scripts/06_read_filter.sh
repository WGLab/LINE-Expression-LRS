#!/bin/bash

# Artifact-Based Filtering Part 1:  Read Filter
# Generate new BAM file of the reads that have 90%  of the read itself within the defined reference L1 Loci and generate a bedgraph of this BAM file

# Set the following parameters in the command:
# 1: Sample Name
# 2: L1 Region Reference File Type 
# 3: Output Files

# Example command:   ./06_read_filter.sh sample_name active 
# (active, inactive, or ORF2)



sample_name=$1

L1_ref_type=$2     # change between: active, inactive, orf2  
echo "Read Filter on the $L1_ref_type Reference L1 Regions"

L1_ref_input=../references/L1Base2_filtered/${L1_ref_type}_filtered.bed
sample_bam_input=../${sample_name}/c_hg38_mapping_LRS/${sample_name}_hg38_mapped.sorted_position.bam


cd ../${sample_name}/d_LINE_quantification/

mkdir -p $L1_ref_type; cd $L1_ref_type

mkdir -p "read_filter"; cd "read_filter"


# Output Files
L1_regions_reads=${sample_name}_L1_regions_reads.bam
echo "Generating filtered BAM file with reads only located within the L1 reference regions..."
samtools view -b -L "$L1_ref_input" -o "$L1_regions_reads" "$sample_bam_input"


read_filter_bam=${sample_name}_read_filter_passed.bam
echo "Removing reads with less than 90% of the read maps to the L1 reference regions"
bedtools intersect -a "$L1_regions_reads" -b "$L1_regions_reads" -f 0.9 > "$read_filter_bam"


sorted_read_filter_bam=${sample_name}_read_filter_passed.sorted_position.bam
echo "Sorting and Indexing the resulting Read Filter BAM file..."
samtools sort "$read_filter_bam" -o "$sorted_read_filter_bam"
samtools index "$sorted_read_filter_bam"







# Generate the bedgraph for the new BAM file that passed the Read Filter 

bedgraph_output=${sample_name}"_bedgraph.bg"         
echo "Generating the bedgraph..." 
bedtools genomecov -ibam "$sorted_read_filter_bam" -bga -split > "$bedgraph_output"


bedgraph_output_clean=${sample_name}"_bedgraph_clean.bg"  
echo "Cleaning the bedgraph..."
grep -v 'fix\|alt\|random\|[(]\|Un' $bedgraph_output > $bedgraph_output_clean

bedgraph_sort_output=${sample_name}"_bedgraph_sorted.bg"  
echo "Sorting the cleaned bedgraph..."
sortBed -i $bedgraph_output_clean > $bedgraph_sort_output

