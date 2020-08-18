# Snakemake pipeline for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

The pipeline is re-implemented and re-wrote by Xin Xu and Xiaoxu Yang, with great help form Martin Breuss, based on a previous version written by Renee D. George, and is maintained by Xin Xu and Xiaoxu Yang.

The pipeline takes a list of paired samples as input, detects mosaic variants and annotates the information needed for downstream analysis and calculate the number of reference counts, alternative counts, mutant allelic fractions (MAFs), and the exact binomial confidence intervals of the MAFs from the "tumor" bam. Details about the input and output are listed below:

----------------------------

## Before start, make sure you have:
#### [Strelka2](https://github.com/Illumina/strelka)
#### [Mutect2 from GATK4](https://github.com/broadinstitute/gatk/releases)
#### [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
#### [BEDTools](https://bedtools.readthedocs.io/en/latest/index.html)
#### [pandas](https://pandas.pydata.org/) and [NumPy](https://numpy.org/) for Python
#### add the [lib](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline/lib) folder to your $PYTHONPATH
#### Be sure you have Python2 on your server.
----------------------------

## Input:
### Below are headers of the input file list
#### sample
User defined name for the specific sample.
#### tumor
The name of the "tumor" sample from the tags in the bam header.
#### normal
The name of the "normal" sample from the tags in the bam header.
#### tumor_path
Path to the "tumor" bam file.
#### normal_path
Path to the "normal" bam file.
#### vcf_path
Path to the variant quality score recalibrated vcfs from haplotype caller. The input should be in vcf format. This will help the calculation such as "near indel".

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
#### strelka_path
Path to your Strelka2 configureStrelkaSomaticWorkflow.py.
#### strelka_config
Path to your Strelka2 configureStrelkaSomaticWorkflow.py.ini.
#### GATK4
Path to your GATK4 jar file.
#### python2
Path to your Python2.
#### annovar
Path to your ANNOVAR annotate_variation.pl.
#### annovar_db
Path to your ANNOVAR databases.
#### java
Path to your java jre.
#### homopolymer_script, ci_script, and summerize_script
Path to the helper scripts, by default the scripts are provided in the helper_scripts folder.
#### bed_file
The total length of the chromosome corresponding to your reference genome file.
#### ref_fasta
Your reference genome.
#### gnomad_af
vcf.gz file with only AF information from gnomAD (corresponding to your reference genome file).
#### ucsc_rpmsk
RepeatMasker track from UCSC genome browser (corresponding to your reference genome file). The default format is #bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEndrepLeft	id
#### repeat_masker
Additional genomic repeats annotated to the output table in bed format (corresponding to your reference genome file).
#### segdup
Segmental duplications annotated to the output table in bed format (corresponding to your reference genome file).

----------------------------

## Output:
### Below are headers of the output table. The output table has the same number of entries as the sum of each of the files from "vcf_path"
#### ID
The name of the "tumor" bam tag where the information of this line is calculated from.
#### CHROM
Chromosome number of each variant.
#### POS
Genomic position of each variant.
#### REF
Reference allele according to the reference genome file given as "ref_fasta".
#### ALT
Alternative allele according to the input variant list.
####  ANNO
Type of variants by ANNOVAR.
#### GENE
Name of gene if the variant is located on/close to the coding region of a gene.
#### GNOMAD_FREQ
Population allele frequency provided in gnomad_af. "0" if not present in the table.
#### REPEAT_MASKER
Repeat region defined in "repeat_masker". "0" if not present in the table.
#### SEGDUP
Segmental duplication defined in "segdup". "0" if not present in the table.
#### HOMOPOLYMER
Homopolymer within ±8 bp sequence. "0" if not present.
####  REF_SEQ
Genomic sequence for variant position ±8 bp sequence (17 bp in total).
#### DINUCLEOTIDE
Dinucleotide repeat within ±8 bp sequence. "0" if not present.
#### NEAR_INDEL
Whether the candidate variant is close to a known germline indel variant by HaplotypeCaller. "0" if not present.
#### UCSC_RPMSK
Annotations in the UCSC repeat masker defined by "ucsc_rpmsk". "pass" if not present.
#### REF_COUNT
Number of reads supporting the reference allele calculated in the "tumor" bam.
#### ALT_COUNT
Number of reads supporting the alternative allele provided in the input vcf calculated in the "tumor" bam.
#### MAF
Mutant allelic fraction calculation from an exact binomial estimation from the "tumor" bam, assuming "REF_COUNT" success form "REF_COUNT"+"ALT_COUNT" trials.
#### LOWER_CI
Lower bound of 95% confidence interval of the exact binomial estimation (Clopper–Pearson interval), assuming "REF_COUNT" success form "REF_COUNT"+"ALT_COUNT" trials.
#### UPPER_CI
Upper bound of 95% confidence interval of the exact binomial estimation (Clopper–Pearson interval), assuming "REF_COUNT" success form "REF_COUNT"+"ALT_COUNT" trials.
#### CI_IS_GREATER
Single char tags for comparison of the variant's "LOWER_CI" with each of the "UPPER_CI" ±5bp of the candidate. "F" if the variant's "LOWER_CI" <= the neighbor's "UPPER_CI" in the "tumor" bam. "P" if the variant's "LOWER_CI" > the neighbor's "UPPER_CI" in the "tumor" bam.
#### NORMAL_REF_COUNT
Number of reads supporting the reference allele calculated in the "normal" bam.
#### NORMAL_ALT_COUNT
Number of reads supporting the alternative allele provided in the input vcf calculated in the "normal" bam.
#### NORMAL_MAF
Mutant allelic fraction calculation from an exact binomial estimation from the "normal" bam, assuming "NORMAL_REF_COUNT" success form "NORMAL_REF_COUNT"+"NORMAL_ALT_COUNT" trials.
#### NORMAL_LOWER_CI
Lower bound of 95% confidence interval of the exact binomial estimation (Clopper–Pearson interval) from the "normal" bam, assuming "NORMAL_REF_COUNT" success form "NORMAL_REF_COUNT"+"NORMAL_ALT_COUNT" trials.
#### NORMAL_UPPER_CI
Upper bound of 95% confidence interval of the exact binomial estimation (Clopper–Pearson interval) from the "normal" bam, assuming "NORMAL_REF_COUNT" success form "NORMAL_REF_COUNT"+"NORMAL_ALT_COUNT" trials.
#### NORMAL_CI_IS_GREATER
Single char tags for comparison of the variant's "NORMAL_LOWER_CI" with each of the "NORMAL_UPPER_CI" ±5bp of the candidate. "F" if the variant's "NORMAL_LOWER_CI" <= the neighbor's "NORMAL_UPPER_CI" in the "normal" bam. "P" if the variant's "NORMAL_LOWER_CI" > the neighbor's "NORMAL_UPPER_CI" in the "normal" bam.
#### TUMOR_IS_BLOOD
Test if the "tumor" sample name given in "tumor" contains the expected key words. Could be edited in the helper script.
#### TUMOR_IS_SPERM
Test if the "normal" sample name given in "normal" contains the expected key words. Could be edited in the helper script.
