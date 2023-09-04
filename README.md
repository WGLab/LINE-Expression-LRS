# LINE Expression Analysis in Long-Read Sequencing 
Bioinformatics pipeline for the identification and quantification of active, full-length Long INterspersed Elements (LINEs)

## Background Information
Long INterspersed Elements (LINEs) are a class of retrotransposons of transposable elements. The abundant number of LINE sequences in the human genome can create genomic structural variants in human populations and lead to somatic alterations in cancer genomes. Evaluating LINE expression will help to better understand the roles they have across diverse tissues and cell lines. With the advent of long-read sequencing technologies, these platforms offer a promising solution by enabling more accurate characterization of LINEs expression profiles. This repository hosts a collection of scripts designed to quantify the expression of LINEs, particularly L1s, from long-read sequencing data. The scripts facilitate data preprocessing, LINE identification, artifact-based filtering, and quantification of L1 expression. 


## Program Requirements
```
1. Minimap2 (v2.24)
2. RepeatMasker (version 4.1.5 released March 2023)
3. Python (v3.7.11)
4. Seqtk (1.2-r94)
5. LongReadSum (v1.3.0)
6. Samtools (1.15.1)
7. Bedtools (v2.30.0)
```
Newer versions for each program may also work. 

## Reference Requirements 
1. [Reference Genome: hg38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.15/)
2. [Comprehensive gene annotation: gencode.v40.annotation.gtf](https://www.gencodegenes.org/human/) (and/or the most recent version)
3. Custom LINE Reference Library: UCSC hg38 RepeatMasker track and LINE/L1 consensus sequences from the Repbase RepeatMasker library (provided in the references folder) 
4. [L1Base2 Active, Inactive, and Intact Only in ORF2 Reference Regions](https://l1base.charite.de/)



# Script Descriptions 
```
01_preprocess_input.sh
02_preprocess_mapping.sh
03_L1_detection.sh
04_map_hg38.sh
05_map_qc_LRS.sh
06_read_filter.sh
07_exon_filter.sh
08_L1_loci_filter.sh
09_map_qc_LRS.sh
10_normalization_wgt_avg.sh
```


## Input File Types
The pipeline is designed to take in a single input file, which can be in either `FASTQ` or `FASTA` format, suitable for both RNA and cDNA sequencing data. If dealing with multiple FASTQ files, it is best to merge the files using the `cat` command beforehand.

If the input is in FAST5 format, it is recommended to use Guppy6 (Available online: https://community.nanoporetech.com) to perform basecalling first, with the default basecalling read Q-score filtering threshold of 7, and then convert the FASTQ files for downstream analysis. 


#### Note: Each example command may take some additional time ran as is, for best results, optimize with the HPC (High-Performance Computing) cluster/computer. 


## 1. Data Preprocessing  Part 1: Initialization 
### `01_preprocess_input.sh`

This script first initializes the subsequent directories and input files required for analysis and was filtered to retain reads longer than 1kb. 

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input  and Reference Files (entire path)
3. Input File Format (FASTQ or FASTA)
4. Data Type (RNA or cDNA)

#### Example command if your input file is in the FASTA format and is from RNA sequencing:   

`./01_preprocess_input.sh sample_name /path/to/input/fasta/file.fasta FASTA RNA`

#### Expected Output in a_dataset folder: 
1. Folder of your specified sample name (specified by parameter 1, no spaces)
2. 4 Sub-folders for downstream analysis
````
a_dataset
b_repeat_masker_process
c_hg38_mapping_LRS
d_LINE_quantification
````
3. Depending on the input parameters set, the script produces the following:

* If the input is in FASTQ format (specified by parameter 3), it converts the input to FASTA format and filters out reads shorter than 1kb.

* If the input data is specified as RNA (parameter 4), it further converts it to cDNA format. If it's already in DNA format, no conversion is performed.

4. FASTA File of your input FASTQ
5. FASTA File with Reads < 1kb Removed





## 2. Data Preprocessing Part 2: Mapping to Custom LINE Reference Library 

### `02_preprocess_mapping.sh`

Next, using a custom LINE reference library as a reference, mapping with Minimap2 was completed as a filtering step to retain reads that contained L1 elements. The custom LINE reference library was constructed by merging the UCSC hg38 RepeatMasker track with the LINE/L1 consensus sequences from the Repbase RepeatMasker library. 

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)

#### Example command: 

`./02_preprocess_mapping.sh sample_name`

#### Expected Output in a_dataset folder:

1. FASTA (.fa) file of the sequences mapped to the custom LINE reference library
2. SAM (.sam) file equivalent 






## 3. L1 Detection Part 1: RepeatMasker

### `03_L1_detection.sh` 

This script first runs RepeatMasker to identify reads with any type of transposable element. The RepeatMasker output was processed to retain reads with less than 10% diverged LINE/L1 elements from the RepeatMasker L1 consensus. 

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)


#### Example command: 

`./03_L1_detection.sh sample_name`

#### Expected Output in b_repeat_masker_process folder:

1. Files specific to the RepeatMasker output: [RepeatMasker Explaination](https://www.repeatmasker.org/webrepeatmaskerhelp.html)

* .out.xm
* .cat
* .align
* .out
* .tbl
* .masked
* .ori.out

2. Using the 

* div10_LINEs.out: Instances from the RepeatMasker output (.out) file with LINE/L1 as the repeat class/family 
* div10_readIDs.txt: Read IDs of the read sequences with LINE/L1 as the repeat class/family
* UHR_ONT_test_div10.fa: FASTA file of these Read IDs 




## 4. L1 Detection Part 2: Mapping Reads with LINE/L1 Elements Detected to the hg38 Reference Genome 

### `04_map_hg38.sh` 

Reads with less than 10% diverged LINE/L1 elements were retained for Minimap2 splice-aware mapping to the hg38 reference genome. 

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)


#### Example command: 

`./04_map_hg38.sh sample_name`

#### Expected Output in c_hg38_mapping_LRS folder:

1. Mapped SAM file to the hg38 reference genome (.sam)
2. Mapped BAM file to the hg38 reference genome (.bam)
3. Sorted by position mapped BAM file to the hg38 reference genome (.bam)
4. Indexed sorted by position mapped BAM file to the hg38 reference genome (.bai)
5. Log output file (.log)



## 5. L1 Detection Part 3: Mapping Statistics with LongReadSum 

### `05_map_qc_LRS.sh` 

LongReadSum was utilized to report mapping statistics for quality check analysis. 

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)








## 6. Artifact-Based Filtering Part 1: Read Filter (Removed)
Artifact-Based Filtering: First, we removed reads where less than 90% of the read itself mapped to the L1 region of interest. With the remaining reads, a bedgraph was generated for downstream calculations of L1 expression levels. 

### `06_read_filter.sh`

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. Specify the Reference L1 Region File (Active, Inactive, or Intact Only in ORF2) 
3. **Full** Path to the Input and Reference Files (entire path)

Example: `06_read_filter.sh sample_name active`










## 7. Artifact-Based Filtering Part 2: L1 Loci with Exon Overlap (Removed)
Next, L1 regions with exon overlap using GENCODE.v43 were removed. The exon annotations from GENCODE.v43 were combined into a new annotation file and using `bedtools intersect`, L1 loci with overlap were removed. 

### `07_exon_filter.sh`

#### To customize the script, provide the following command-line parameters:
1. **Full** Path to the Input and Reference Files (entire path)









## 8. Artifact-Based Filtering Part 3: L1 Loci Filter (Removed) 
Next, regions with fewer than two mapped reads, inconsistent read start position among the reads located within the L1 region (threshold 100bps), and overall starting positions too far from the L1 region start position (threshold 1.5kb) were removed from analysis to filter out false positives.  

### `08_L1_loci_filter.sh`

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. **Full** Path to the Input and Reference Files (entire path)

The final number of reads for each region is located in the file titled: `active_coverage_for_weighted_avg.bed` (for example, if you are calculating over the active L1 regions)








## 9. Mapping Statistics with LongReadSum on Filtered BAM File

### `09_map_qc_LRS.sh` 

LongReadSum was utilized to report mapping statistics for quality check analysis. 

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. Type of L1 Category you are analyzing (active, inactive, ORF2) 
3. **Full** Path to the Input and Reference Files (entire path)









## 10. L1 Quantification: Normalization by Total Number of Reads and Weighted Average
To determine the general expression level of these L1 loci, the calculated coverage values were first normalized by the total number of reads in the sample. Then, a weighted average was computed across all L1 loci in each of the reference L1Base2 categories: active, inactive, and intact only in ORF2, for each sample. This weighted average involved assigning weights based on the length of the reference L1 element with its corresponding expression value, scaled in millions, as displayed below. For i takes the L1 reference category, "Active," "Inactive," and "Intact only in ORF2" to represent the different categories. 
<img width="1179" alt="Screenshot 2023-09-01 at 5 08 03 PM" src="https://github.com/WGLab/LINE-Expression-LRS/assets/89222332/6582f5e5-105b-462a-80d6-0164fbae5bf6">

### `10_normalization_wgt_avg.sh`

#### To customize the script, provide the following command-line parameters:
1. Sample Name  (all one word)
2. Type of L1 Category you are analyzing (active, inactive, ORF2) 
3. **Full** Path to the Input and Reference Files (entire path)








## Contact
If you have any questions/issues/bugs, please post them on GitHub. 

## Reference
Rybacki K, Xia M, Ahsan MU, Xing J*, Wang K*. Assessing the expression of Long INterspersed Elements (LINEs) via long-read sequencing in diverse human tissues and cell lines. *Submitted*, 2023

