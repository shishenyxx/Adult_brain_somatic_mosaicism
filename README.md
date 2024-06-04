## Adult_brain_somatic_mosaicism

This repository is a collaboration of codes and scripts used in the analysis of [a somatic mosaicism study of the NIMH Brain Somatic Mosaicism Network](https://bsmn.synapse.org/Explore/Studies/DetailsPage?id=syn22269661) on the entire human cortical regions, subregions, sorted populations, and single nuclei from ID01(Raw data are available at [the NDA website under study #919](https://nda.nih.gov/study.html?id=919)) and ID02-04 (Raw data are available [here](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA736951&o=acc_s%3Aa)). The 300x WGS panel of normal data is also available [here](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA660493&o=acc_s%3Aa).

-----------------------------------

### 1. Pipelines for the process of whole-genome sequencing data
#### 1.1 Pipelines for WGS data pre-process, alignment, post-process, and quality control
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_processing_pipeline) for converting the Dragen bam into fastqs, mapping them with new parameters to GRCh37, and finish the indel re-alignment as well as base quality score recalibration.

[Codes and scripts](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/QC) for WGS quality control.

#### 1.2 Pipelines for mosaic SNV/indel calling and variant annotations
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

Pipelines for [MuTect2 (single mode)](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode), followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

PBS script for [MosaicHunter (single mode)](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicHunter_single_mode_pipeline), followed by the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

After variant calling from different strategies, variants were annotated and filtered by [a python script](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/pipelines/WGS_SNV_indel_calling_pipeline/Filter_and_annotate_candidate_mosaic_variants.py) and positive mosaic variants as well as the corresponding tissue and additional information were annotated.

-----------------------------------

### 2. Pipelines for the process of Massive Parallel Amplicon Sequencing (MPAS) and single-nuclei MPAS (snMPAS)
#### 2.1 Pipelines for MPAS data alignment and processing
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/MPAS_and_snMPAS_processing_pipeline) for alignment, processing, and germline variant calling of MPAS and snMPAS reads.

#### 2.2 Pipelines for AF quantification and variant annotations
[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

-----------------------------------

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting for ID01
[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/pipelines/Codes_for_mosaic_variant_annotations_after_MPAS.py) to filter and annotate on MPAS and snMPAS data.

[Codes and config files](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/circos) for the Circos plot of square root-transformed AFs measured by MPAS.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/permutation) for permutation analysis from gnomAD and [codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Genomic_enrichment/Plot_enrichment.r) for plotting the permutation result.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Supplement_data_single_cell/AF_for_single_cell.r) for plotting of AF measured in snMAPS.

#### 3.2 Pipelines for mosaic variant determination, annotations, and plotting for ID02, 03, and 04
[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/4dbsm_final_data_vali_new_data11.py) to filter and annotate and plot for ID02, 03, and 04.

#### 3.3 Pipelines for statistically analysis, QC, and the related plotting
Codes for the QC of [MPAS](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/QC/Plot_MPAS_het_and_ref_homo_controls.r) and [snMPAS](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/QC/Plot_snMPAS_het_and_ref_homo_controls.r) based on the heterozygous and reference homozygous control variants in the panel.

Codes for sorted population, and the QC for sorting were already described in the [previous publication](https://science.sciencemag.org/content/366/6469/1134/tab-pdf).

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Correlation_clustering/corr_clustermap.py) for correlation analysis and cluster representation of the AFs measured by MPAS.

[Codes and permutations](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/sub_lobar_permutation) for the analysis for the distribution of sublobar areas.

[Computational simulations and plotting](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/plotting/Left_right_founder_estimation) for left-right starting populations.

### 4. [Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/Codes_for_plotting_main_figure_panels_based_on_MPAS_and_snMPAS_annotation.py) for the plotting of panels in the main figures and supplements for ID01

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/plotting/4dbsm_final_data_vali_new_data_LIBD02.py) for plotting supplement panels for ID02.

### 5. [Jupyter Notebook](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/umap/4dbsm_umap.ipynb) for the UMAP
UMAP clustering from MPAS from bulk and sorted samples from ID01.

### 6. [Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/Lineage_construction/lineage_construction.py) for the lineage construction
Lineage reconstruction based on snMPAS and MPAS evaluated genotypes from sorted nuclei as well as bulk and sorted populations from ID01.

### 7. Contact information

:email: Xiaoxu Yang: [u6055394@utah.edu](mailto:u6055394@utah.edu), [xiaoxuyanglab@gmail.com](mailto:xiaoxuyanglab@gmail.com), [xiy010@health.ucsd.edu](mailto:xiy010@health.ucsd.edu)

:email: Martin Breuss: [martin.breuss@cuanschutz.edu](mailto:martin.breuss@cuanschutz.edu)

:email: Joseph Gleeson: [jogleeson@health.ucsd.edu](mailto:jogleeson@health.ucsd.edu)

### 8. Cite the codes

<b>Breuss MW, Yang X, Antaki D, Schlachetzki JCM, <i>et al.,</i> Gleeson JG. [Somatic mosaicism reveals clonal distributions of neocortical development.](https://www.nature.com/articles/s41586-022-04602-7) 2022. (<i>Nature</i>, DOI:[10.1038/s41586-022-04602-7](https://doi.org/10.1038/s41586-022-04602-7),  PMID:[35444276](https://pubmed.ncbi.nlm.nih.gov/35444276/))</b>


<img src="https://user-images.githubusercontent.com/17311837/223878022-1b92c9b1-6a9b-481e-be1d-6d1fbdd66668.png" alt="Sperm_Mosaic_Cover" width=40%> 

