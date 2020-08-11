## Adult_brain_somatic_mosaicism

This repository is a collaboration of codes and scripts used in the analysis of a somatic mosaicism [study](https://staging.bsmn.synapse.org/Explore/Studies/DetailsPage?id=syn22269661) on an entire human coretical regions, subregions, sorted populations and single nuclei. Raw data is available at the NDA website under study [#919](https://nda.nih.gov/study.html?id=919).

Cite the codes: [Breuss, Yang, Antaki, Schlachetzki et al. Somatic mosaicism in the mature brain reveals clonal cellular distributions during cortical development. 2020. BioRxiv.](https://doi.org/10.1101/2020.08.10.244814)

-----------------------------------

### 1. Pipelines for the process of whole-genome sequencing data
#### 1.1 Pipelines for WGS data pre-process, alignment, post-process, and quality control
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_processing_pipeline) for converting the Dragen bam into fastqs, mapping them with new parameters to GRCh37, and finish the indel re-alignment as well as base quality score recalibration.

[Codes and scripts](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/QC) for WGS quality control.

#### 1.2 Pipelines for mosaic SNV/indel calling and variant annotations
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

Pipelines for [MuTect2 single mode](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode), followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

PBS script for [MosaicHunter single mode](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicHunter_single_mode_pipeline), followed by the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

After variant calling from different strategies, variants were annotated and filtered by [a python script](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/pipelines/WGS_SNV_indel_calling_pipeline/Filter_and_annotate_candidate_mosaic_variants.py) and positive mosaic variants as well as the corresponding tissue and additional information were annotated.

-----------------------------------

### 2. Pipelines for the process of Massive Parallel Amplicon Sequencing (MPAS) and single-nuclei MPAS (snMPAS)
#### 2.1 Pipelines for MPAS data alignment and processing
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/MPAS_and_snMPAS_processing_pipeline) for alignment, processing, and germline variant calling of MPAS and snMPAS reads.

#### 2.2 Pipelines for AF quantification and variant annotations
[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

-----------------------------------

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting
[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/pipelines/Codes_for_mosaic_variant_annotations_after_MPAS.py) to filter and annotate on MPAS and snMPAS data.

[Codes and config files](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/circos) for the Circos plot of square root-transformed AFs measured by MPAS.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/permutation) for permutation analysis from gnomAD and [codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Genomic_enrichment/Plot_enrichment.r) for plotting the permutation result.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Supplement_data_single_cell/AF_for_single_cell.r) for plotting of AF measured in snMAPS.
#### 3.2 Pipelines for statistically analysis, QC, and the related plotting
Codes for the QC of [MPAS](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/QC/Plot_MPAS_het_and_ref_homo_controls.r) and [snMPAS](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/QC/Plot_snMPAS_het_and_ref_homo_controls.r) based on the heterozygous and reference homozygous control variants in the panel.

Codes for sorted population, and the QC for sorting were already described in the [previous publication](https://science.sciencemag.org/content/366/6469/1134/tab-pdf).

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Correlation_clustering/corr_clustermap.py) for correlation analysis and cluster representation of the AFs measured by MPAS.

[Computational simulations and plotting](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/Left_right_founder_estimation) for left-right starting populations.

## 4 [Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Codes_for_plotting_main_figure_panels_based_on_MPAS_and_snMPAS_annotation.py) for the plotting of panels in the main figures and supplements

## 5 [Jupyter Notebook](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/umap/4dbsm_umap.ipynb) for the UMAP

## 6 [Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/Lineage_construction/lineage_construction.py) for the lineage construction
