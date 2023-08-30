# LINE Expression Analysis in Long-Read Sequencing 
Bioinformatics pipeline for the identification and quantification of active, full-length Long INterspersed Elements (LINEs)

## Background Information

Long INterspersed Elements (LINEs) are a class of retrotransposons of transposable elements. The abundant number of LINE sequences in the human genome can create genomic structural variants in human populations and lead to somatic alterations in cancer genomes. Evaluating LINE expression will help to better understand the roles they have across diverse tissues and cell lines. With the advent of long-read sequencing technologies, these platforms offer a promising solution by enabling more accurate characterization of LINEs expression profiles. This repository hosts a collection of scripts designed to quantify the expression of LINEs, particularly L1s, from long-read sequencing data. The scripts facilitate data preprocessing, LINE identification, artifact-based filtering, and quantification of L1 expression. 


## Program Requirements
```
Minimap2
RepeatMasker
Python
Seqtk
LongReadSum
Samtools
Bedtools
```

# Script Descriptions 

## Input File Types
The pipeline is designed to take in a single input file, which can be in either `FASTQ` or `FASTA` format, suitable for both RNA and cDNA sequencing data. If dealing with multiple FASTQ files, it is best to merge the files using the `cat` command beforehand.


## 1. Data Preprocessing  Part 1: Initialization 

### `a_preprocess.sh`

This script first initializes the subsequent directories and input files required for analysis and was filtered to retain reads longer than 1kb. 

To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input  and Reference Files (entire path)
3. Input File Format (FASTQ or FASTA)
4. Data Type (RNA or cDNA)


## 2. Data Preprocessing Part 2: Mapping to Custom LINE Reference Library 

### `a_preprocess_mapping.sh`

Next, using a custom LINE reference library as a reference, mapping with Minimap2 was completed as a filtering step to retain reads that contained L1 elements. The custom LINE reference library was constructed by merging the UCSC hg38 RepeatMasker track with the LINE/L1 consensus sequences from the Repbase RepeatMasker library. 

To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)


## 3. L1 Detection Part 1: RepeatMasker

### `b_L1_detection.sh` 

This script first runs RepeatMasker to identify reads with any type of transposable element. The RepeatMasker output was processed to retain reads with less than 10% diverged LINE/L1 elements from the RepeatMasker L1 consensus. 

To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)


## 4. L1 Detection Part 2: Mapping Reads with LINE/L1 Elements Detected to the hg38 Reference Genome 

### `b_map_hg38.sh` 

Reads with less than 10% diverged LINE/L1 elements were retained for Minimap2 splice-aware mapping to the hg38 reference genome. 

To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)


## 5. L1 Detection Part 3: Mapping Statistics with LongReadSum 

### `b_map_qc_LRS.sh` 

LongReadSum was utilized to report mapping statistics for quality check analysis. 

To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)


## 6. Artifact-Based Filtering
Artifact-Based Filtering: First, we removed reads where less than 90% of the read itself mapped to the L1 region of interest. Then, L1 regions with exon overlap, fewer than two mapped reads, inconsistent read start position among the reads within the L1 region (threshold 100bps), and overall starting positions too far from the L1 region start position (threshold 1.5kb) were removed from analysis to filter out false positives. 




## 5. L1 Quantification





## Contact
If you have any questions/issues/bugs, please post them on GitHub. 

## Reference
