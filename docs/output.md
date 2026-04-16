# Output

## Directory structure

```
results/
├── pipeline_info/
│   ├── validated_samplesheet.csv
│   ├── validated_contrasts.csv
│   └── software_versions.yml
├── qc/
│   ├── qc_metrics.csv
│   ├── qc_vst_matrix.tsv
│   ├── qc_library_sizes.{pdf,png}
│   ├── qc_count_distribution.{pdf,png}
│   ├── qc_pca.{pdf,png}
│   ├── qc_correlation_heatmap.{pdf,png}
│   └── qc_dispersion.{pdf,png}
├── normalization/
│   ├── deseq2_dds.rds
│   ├── deseq2_vst_counts.tsv
│   └── deseq2_size_factors.csv
├── dge/
│   └── <contrast_id>/
│       ├── <contrast_id>_DESeq2_results.csv
│       ├── <contrast_id>_volcano.{pdf,png}
│       ├── <contrast_id>_MA.{pdf,png}
│       └── <contrast_id>_top_genes_heatmap.{pdf,png}
├── gsea/
│   └── <contrast_id>/
│       ├── <contrast_id>_GO_BP_gsea.csv
│       ├── <contrast_id>_GO_MF_gsea.csv
│       ├── <contrast_id>_GO_CC_gsea.csv
│       ├── <contrast_id>_GO_all_gsea.csv
│       ├── <contrast_id>_KEGG_gsea.csv
│       ├── <contrast_id>_Hallmarks_gsea.csv
│       ├── <contrast_id>_Reactome_gsea.csv
│       └── *_dotplot.{pdf,png}
└── dashboard/
    └── (Shiny app source + staged data)
```

## Key files

### QC
- `qc_metrics.csv` — library size, detected genes, size factors per sample.
- `qc_vst_matrix.tsv` — blind VST-transformed counts (for exploratory visualisation).
- `qc_pca.png` — PCA of all samples coloured by condition.
- `qc_correlation_heatmap.png` — sample-sample Pearson correlation.

### Normalization
- `deseq2_dds.rds` — fitted `DESeqDataSet` object (loaded by DGE processes).
- `deseq2_vst_counts.tsv` — model-aware VST matrix (non-blind) for all samples.
- `deseq2_size_factors.csv` — DESeq2 size factors.

### DGE
- `*_DESeq2_results.csv` — full results table with columns: `gene_id`, `gene_name`, `baseMean`, `log2FoldChange` (ashr-shrunken), `lfcSE`, `stat`, `pvalue`, `padj`.
- `*_volcano.png` — volcano plot with top 20 significant genes labelled.
- `*_MA.png` — MA plot.
- `*_top_genes_heatmap.png` — heatmap of top 50 genes (by padj).

### GSEA
- Each CSV contains one row per gene set with: `Description`, `NES`, `pval`, `padj`, `size`, `core_enrichment`/`leadingEdge`.
- `*_GO_all_gsea.csv` — combined GO results (BP + MF + CC) with `ontology` column.
- Dot plots show top 20 significant gene sets coloured by adj. p-value.

### Dashboard
Run `bash dashboard/launch_dashboard.sh results/` to launch the interactive Shiny app locally.
