# Snakemake pipeline for MosaicForecast mosaic variant classifier from WGS data, based on the result of MuTect2 single mode

The pipeline is re-implemented from the official MosaicForecast pipeline by Yanmei Dou.

The pipeline takes a list of bam samples and vcf files from MuTect2 single mode as input, detects mosaic variants for downstream analysis.
----------------------------

## Before start, make sure you have:
#### [MosaicForecast](https://github.com/parklab/MosaicForecast)
#### [pandas](https://pandas.pydata.org/) and [NumPy](https://numpy.org/) for Python

----------------------------

## Input:
### Below are headers of the input file list
#### sample
User defined name for the specific sample.
#### tumor_path
Path to the input bam file.
#### Mutect_SM_vcf_file
Path to the vcf.gz file generated from MuTect2 single mode.

----------------------------

## Config files:
### Below are files you need to prepare for the annotation scripts, saved in the file snake_conf.yaml
#### input_files
Path to the input file list.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqal to two- or three-times the number of total listed variants.
#### ref_fasta
Your reference genome.
#### Other parameters and output format please refer to [MosaicForecast](https://github.com/parklab/MosaicForecast).

