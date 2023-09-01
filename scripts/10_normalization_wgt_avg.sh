#!/bin/bash

# Normalization by Total Number of Reads from LongReadSum, followed by a Weighted Average

# Set the following parameters in the command:
# 1: Sample Name
# 2: L1 Region Reference File Type
# 3: Output Files

# Example command:   ./10_normalization_wgt_avg.sh sample_name active
# (active, inactive, or ORF2)


sample_name=$1
L1_ref_type=$2     # change between: active, inactive, orf2


cd ../${sample_name}/d_LINE_quantification/${L1_ref_type}/

echo "Normalization by Total Number of Reads..."
qc_report_input=../${sample_name}/d_LINE_quantification/${L1_ref_type}/read_filter/${sample_name}/bam_summary.txt

# function to calculate coverage

# Extract the number of read value from summary.txt
value=$(awk -F '\t' 'NR==1 {print $2}' "$qc_report_input")

#value=$(awk '{print $5}' $qc_report_input)
echo "Number of reads: $value"


input_ref_cov=../${sample_name}/d_LINE_quantification/${L1_ref_type}/L1_loci_filter/${L1_ref_type}_coverage_for_weighted_avg.bed




# output files
normalized_ref_regions="normalized_${L1_ref_type}_regions.bed"

awk -v divisor="$value" -v OFS="\t" '{$NF = (divisor != 0) ? $NF / divisor : 0; print}' "$input_ref_cov" > "$normalized_ref_regions"

echo "Normalization Complete!"




# Part 2: Weighted Average Calculation 
echo "Weighted Average Calculation..." 

# output file 
weighted_average_cov="coverage_weighted_avg.bed"


weighted_sum=0
total_weight=0

while IFS=$'\t' read -r line || [[ -n "$line" ]]; do
    elements=($line)
    start=${elements[1]}
    end=${elements[2]}
    coverage=$(echo "${elements[-1]}" | awk '{print $NF}')
    region_size=$((end - start + 1))
    weighted_sum=$(awk "BEGIN {print $weighted_sum + ($coverage * $region_size)}")
    total_weight=$((total_weight + region_size))
done < "$input_ref_cov"

if [ "$total_weight" -ne 0 ]; then
    weighted_average=$(awk "BEGIN {print $weighted_sum / $total_weight}")
else
    weighted_average=0
fi

# Write the Sample Name and Coverage Values to the new file 
echo -e "Sample Name\tWeighted Average" > "$weighted_average_cov"
echo -e "${sample_name}\t$weighted_average" >> "$weighted_average_cov"

echo "Calculations for ${sample_name}, over the ${L1_ref_type} regions is complete!"




