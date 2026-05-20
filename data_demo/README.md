# bulkAnnex Demo Dataset

Synthetic bulk RNA-seq dataset mimicking nf-core/rnaseq output (star_salmon pipeline).
Designed to exercise the full bulkAnnex pipeline without using patient or sensitive data.

## Dataset

| Property | Value |
|---|---|
| Samples | 12 (2 norm_groups × 3 Scramble + 3 KD) |
| norm_groups | CellLine_A, CellLine_B |
| Genes | 1,998 human (real Ensembl IDs + HGNC symbols) |
| Biological scenario | RNA helicase knockdown (DDX/DHX family) |
| KD targets | DDX3X, DDX5, DDX17, DDX21, DHX9, DHX36 (3–5× down) |
| Secondary effects | Apoptosis genes up (~2×), cell cycle genes down (~2×) |
| Expected GSEA pathways | Cell cycle, apoptosis, RNA processing, PI3K-AKT signaling |
| CellLine_B difference | ~30% stronger KD effect, ~600 genes shifted at baseline |

## Setup

**Step 1 — Generate the counts matrix** (run once):

```bash
python data_demo/generate_demo.py
```

This creates:
```
data_demo/salmon.merged.gene_counts.tsv
```

**Step 2 — Run the pipeline:**

```bash
# Local (conda environment)
bash data_demo/run_command.sh

# HPC / SLURM (Apocrita) — run from the project root
sbatch data_demo/submit.sh
```

## Expected results

`results_demo/` will contain:

- **QC plots** — library sizes, count distribution, PCA (CellLine_A and CellLine_B cluster separately), correlation heatmap
- **Normalization** — VST-normalized counts, size factors per sample
- **DGE** — two contrast directories:
  - `CellLine_A_KD_vs_Scramble/` — volcano plot with DDX/DHX genes clearly down
  - `CellLine_B_KD_vs_Scramble/` — same program, stronger fold-changes
- **GSEA** — enrichment in RNA helicase / RNA processing, apoptosis, and cell cycle pathways
- **Dashboard** — interactive exploration of all results
