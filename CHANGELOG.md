# Changelog

All notable changes to bulkAnnex will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- QC PCA is now also produced per `norm_group` (`qc_pca_<norm_group>.{pdf,png}`) alongside the global `qc_pca_all.{pdf,png}`.
- Dashboard QC > PCA: default colouring by sample name, a `norm_group` selector, and PNG/PDF download buttons.
- Initial pipeline implementation (v1.0.0-dev)
- INPUT_CHECK module: samplesheet and contrasts validation
- DESEQ2_QC module: library sizes, PCA, correlation heatmap, dispersion
- DESEQ2_NORMALIZATION module: VST normalization, filtered counts, fitted DDS
- DESEQ2_DGE module: per-contrast differential expression with volcano, MA, heatmap plots
- GSEA_ANALYSIS module: GO (BP/MF/CC), KEGG, Reactome enrichment; ridgeplots; pathview KEGG maps
- SHINY_DASHBOARD module: interactive exploration of all results
- Support for human (GRCh38) and mouse (GRCm39)
- Auto-detection of batch variable in DESeq2 design formula
- Profiles: singularity, docker, conda, test, slurm, apocrita
- Demo dataset: 15 NB4 samples (5 conditions × 3 experiments)

### Changed
- QC PCA plots are now coloured by sample name (legend = sample) instead of condition; `qc_pca` was renamed to `qc_pca_all`.
