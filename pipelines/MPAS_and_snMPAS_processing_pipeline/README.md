# Snakemake pipeline for alignment and germline heterozygous variant calling from MPAS and snMPAS data

The pipeline is maintained by Chen li, Xiaoxu Yang, and Renee George 

----------------------------

## Before start, make sure you have:
#### [BWA](http://bio-bwa.sourceforge.net/) for mapping.
#### [SAMtools](http://www.htslib.org/) for post-alignment processing.
#### [BCFtools](http://samtools.github.io/bcftools/bcftools.html) now part of SAMtools.
#### [java](https://www.java.com/en/download/help/linux_x64_install.xml) for Linux.
#### [Picard Tools](https://broadinstitute.github.io/picard/) from BROAD institute.
#### [GATK](https://github.com/broadgsa/gatk/releases) 3.8 for this pipeline.

----------------------------

## Input:
### Below are headers of the input file sample_list.txt
#### NAME
Unique IDs for any input Dragen sample read pair.
#### R1.fastq.gz
Read1 fastq sequences from the sequencer.
#### R2.fastq.gz
Read2 fastq sequences from the sequencer.

----------------------------

## Config files:
### Configs of this pipeline cound be found in config.yaml

----------------------------

## Output:
The final output bam is in the recaled_bams folder.

