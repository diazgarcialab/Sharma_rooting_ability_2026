# Sharma_rooting_ability_2026

Code and data for:

> Sharma S, Lupo Y, Munoz J, Cochetel N, Nunez V, Gaspar A, Torres-Lomas E, Cantu D,
> Diaz-Garcia L. *Harnessing Vitis diversity to dissect and predict adventitious rooting
> traits in grapevine.* The Plant Genome (under review, TPG-2026-08-0296).

Corresponding author: Luis Diaz-Garcia, diazgarcia@ucdavis.edu

This repository is provided as is, to document the analysis behind the paper. It is a
research pipeline, not a software package: it is not maintained, no support is offered,
and it has only been run on the authors' own machines. Everything needed to reproduce
the reported figures, tables and statistics from the deposited data is here, with the
exceptions listed under **What is not included**.


## What is not included, and where to get it

- **GEMMA association output** (`02_clean_data/gwas_results/*.assoc.txt`, about 390 MB
  per trait) exceeds GitHub's file limits. It is regenerated from the VCF with
  `01_code/gemma_rooting.sh`, which documents the PLINK and GEMMA commands and
  filters used.
- **Sequencing reads** are in the NCBI Sequence Read Archive under BioProjects
  PRJNA731597, PRJNA842753 and PRJNA984685. The joint-called VCF is derived from these
  as described in the Methods.
- **Reference genome and annotation** (*V. girdiana* SC2 v1.0, haplotype 1) are
  available from www.grapegenomics.com.
- **WorldClim 2.1 rasters** are downloaded on first run by `02_figure_01_accession_map.R`
  into `00_raw_data/spatial_cache/`.
