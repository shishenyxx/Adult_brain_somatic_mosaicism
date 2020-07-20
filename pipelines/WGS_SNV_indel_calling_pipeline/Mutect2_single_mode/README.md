# Snakemake pipeline for MuTect2 (single mode)

The pipeline is wrote by Xin Xu, and is maintained by Xin Xu and Xiaoxu Yang.

The pipeline takes a list of samples assigned with "tumor" or "normal" (panel of normals) as input, mosaic SNV/indels from each "tumor" bam were called with consideration of the information from the panel of normals. Details about the input and output are listed below:

----------------------------

## Before start, make sure you have:
#### [Mutect2 from GATK4](https://github.com/broadinstitute/gatk/releases)
#### [Java](https://java.com/en/download/help/linux_x64_install.xml)

----------------------------

## Input:
### Below are headers of the input file list
#### ID
User defined name for the specific sample.
#### GROUP
The group of the sample, either "tumor" or "normal".
#### BAM_PATH
Path to the bam file.

----------------------------

## Config files:
### Below are files you need to prepare for the annotation scripts, saved in the file snake_conf.yaml
#### input_files
Path to the input file list.
#### n_intervals
Number of sub jobs to split for each sample. Jobs will be submitted independently to the server.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqal to two- or three-times the number of total listed variants.
#### GATK4
Path to your GATK4 jar file.
#### java
Path to your java jre.
#### bed_file
The total length of the chromosome corresponding to your reference genome file.
#### ref_fasta
Your reference genome.
#### gnomad
vcf.gz file with only AF information from gnomAD (corresponding to your reference genome file).

----------------------------

## Output:
Raw vcfs and filtered vcfs from Mutect2. The raw vcfs will be used in followup methods.
