# Snakemake pipeline for preprocessing of the DRAGEN bam and remapping them to GRCh37 genome for the follow up mosaic variant calling algorithms

The pipeline is originally wirtten by Danny Antaki, and re-implemented by Xin Xu, with help form Martin Breuss and Xiaoxu Yang

----------------------------

## Before start, make sure you have:
#### [java](https://www.java.com/en/download/help/linux_x64_install.xml) for Linux.
#### [Picard Tools](https://broadinstitute.github.io/picard/) from BROAD institute.
#### [BWA](http://bio-bwa.sourceforge.net/) for mapping.
#### [SAMtools](http://www.htslib.org/) for post-alignment processing.
#### [GATK](https://github.com/broadgsa/gatk/releases) 3.8 for this pipeline.
#### [Sambamba](https://lomereiter.github.io/sambamba/) to help with the fast processing.

----------------------------

## Input:
### Below are headers of the input file list
#### UNIQ_ID
Unique IDs for any input Dragen bam file.
#### SAMPLE_ID
The ID of the sample, might be shared between several bams.
#### BAM_PATH
The path to the input Dragen bam.

----------------------------

## Config files:
### Below are files you need to prepare for the annotation scripts, saved in the file snake_conf.yaml
#### input_files
Path to the input file list.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqal to two- or three-times the number of total listed variants.
#### temp_dir
Path to the temporary files. The temporary files will be deleted after the process is successfully finished.

#### java
Path to your java jre.
#### picard
Path to your picard.jar.
#### bwa
Path to your BWA.
#### samtools
Path to your SAMtools.
#### gatk
Path to your GenomeAnalysisTK.jar
#### sambamba
Path to your executable sambamba

#### hg19_fasta
Your reference genome.
#### mills_indel
vmills indel vcf file (corresponding to your reference genome file).
#### gp1000_indel
List for common indels from the 1000 Genome Project in vcf format (corresponding to your reference genome file).
#### db_gap
Snp list of your dbsnp file in vcf format (corresponding to your reference genome file).


----------------------------

## Output:
The final output bam is in the recaled_bams folder.


