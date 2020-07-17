## Adult_brain_somatic_mosaicism

This repository is a collaboration of codes and scripts used in the analysis of human somatic mosaic variant study.

-----------------------------------

### 1. Pipelines for the process of whole-genome sequencing data
#### 1.1 Pipelines for WGS data pre-process, alignment, post-process, and quality control
[Codes and scripts](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/QC) for WGS quality control.
#### 1.2 Pipelines for mosaic SNV/indel calling and variant annotations

-----------------------------------

### 2. Pipelines for the process of Massive Parallel Amplicon Sequencing (MPAS) and single-nuclei MPAS (snMPAS)
#### 2.1 Pipelines for MPAS data alignment and processing
#### 2.2 Pipelines for AF quantification and variant annotations
[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

-----------------------------------

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting
[Codes and config files](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/circos) for the Circos plot of square root-transformed AFs measured by MPAS.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/permutation) for permutation analysis from gnomAD and [codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Genomic_enrichment/Plot_enrichment.r) for plotting the permutation result.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Supplement_data_single_cell/AF_for_single_cell.r) for plotting of AF measured in snMAPS.
#### 3.2 Pipelines for statistically analysis, QC, and the related plotting
Codes for QC for [MPAS](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/QC/Plot_MPAS_het_and_ref_homo_controls.r) and [snMPAS](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/QC/Plot_snMPAS_het_and_ref_homo_controls.r) based on the heterozygous and reference homozygous control variants in the panel.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Correlation_clustering/corr_clustermap.py) for correlation analysis and cluster representation of the AFs measured by MPAS.

[Computational simulations and plotting](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/Left_right_founder_estimation) for left-right starting populations.
#### 3.3 Pipelines for sorted population, the QC of sorting
#### 3.4 Pipelines for further analysis

Umap clustering for variants and sequenced samples
