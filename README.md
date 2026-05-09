# bulkAnnex

**Downstream bulk RNA-seq analysis pipeline** — QC → normalization → differential expression → GSEA → interactive Shiny dashboard.

Takes the gene count matrix from [nf-core/rnaseq](https://nf-co.re/rnaseq/) (`salmon.merged.gene_counts.tsv`) as input.

---

## Quick start

```bash
nextflow run main.nf -profile test,singularity
```

---

## Inputs

| Parameter | Description |
|-----------|-------------|
| `--input` | `samplesheet.csv` — columns: `sample`, `condition`, `batch` (optional), `norm_group` (optional) |
| `--counts` | `salmon.merged.gene_counts.tsv` — gene_id, gene_name, one column per sample |
| `--contrasts` | `contrasts.csv` — columns: `contrast_id`, `variable`, `reference`, `treatment`, `norm_group` (optional) |
| `--organism` | `human` (default) or `mouse` |
| `--outdir` | Output directory (default: `results`) |

Copy-paste templates: [assets/samplesheet_template.csv](assets/samplesheet_template.csv) · [assets/contrasts_template.csv](assets/contrasts_template.csv)  
Full worked example: [data_demo/samplesheet.csv](data_demo/samplesheet.csv) · [data_demo/contrasts.csv](data_demo/contrasts.csv)  
Column reference: [docs/usage.md](docs/usage.md)

---

## Group-aware normalisation (`norm_group`)

When an experiment mixes biologically distinct sample types (e.g. cell lines and primary patient cells), fitting a single DESeq2 model conflates their dispersion estimates. The optional `norm_group` column lets you run independent normalisation per group while keeping all downstream steps (DGE, GSEA, dashboard) unified.

### How it works

Add a `norm_group` column to both `samplesheet.csv` and `contrasts.csv`:

```csv
# samplesheet.csv
sample,condition,norm_group
NB4_Exp4_Scramble,NB4_Scramble,NB4
NB4_Exp4_DDX41sh1,NB4_DDX41_sh1,NB4
MSCline_Exp8_Scramble,MSCline_Scramble,MSCline
HealthyDonor_MNC_H1,MNC_Healthy,MNC
DDX41Patient_MNC_P1,MNC_DDX41_Patient,MNC

# contrasts.csv
contrast_id,variable,reference,treatment,norm_group
NB4_DDX41_sh1_vs_Scramble,condition,NB4_Scramble,NB4_DDX41_sh1,NB4
MSCline_DDX41_sh1_vs_Scramble,condition,MSCline_Scramble,MSCline_DDX41_sh1,MSCline
MNC_DDX41_Patient_vs_Healthy,condition,MNC_Healthy,MNC_DDX41_Patient,MNC
```

- Each `norm_group` runs DESeq2 on its own samples only → dispersion estimated within-group.
- QC (PCA, library sizes, correlation heatmap) always uses **all samples together** for a global overview.
- `dge/` and `gsea/` output directories remain flat — contrast IDs are unique across groups.
- **Backward compatible**: omitting `norm_group` is equivalent to `norm_group = "all"` for every row; output structure is identical to previous runs.

### Output structure with norm_group

```
results/
├── qc/                          # global QC — all samples
├── normalization/
│   ├── NB4/                     # DESeq2 model for NB4 samples only
│   │   ├── deseq2_dds.rds
│   │   └── deseq2_vst_counts.tsv
│   ├── MSCline/
│   ├── MNC/
│   └── MSC/
├── dge/                         # flat — one dir per contrast_id (unchanged)
└── gsea/
```

Without `norm_group` (or all rows set to `"all"`):

```
results/
├── normalization/               # flat — single DESeq2 model (previous behaviour)
│   ├── deseq2_dds.rds
│   └── deseq2_vst_counts.tsv
├── dge/
└── gsea/
```

---

## Pipeline steps

1. **Input validation** — checks samplesheet, contrasts, and counts consistency
2. **QC** — library sizes, count distributions, PCA, correlation heatmap, dispersion
3. **Normalization** — DESeq2 VST, fitted DDS saved for DGE reuse
4. **DGE** — DESeq2 results per contrast, volcano, MA, top-gene heatmap
5. **GSEA** — GO (BP/MF/CC), KEGG, Reactome via clusterProfiler/ReactomePA + pathview
6. **Dashboard** — self-contained Shiny app with 5 interactive tabs

---

## Usage

```bash
# HPC (Apocrita/SLURM) — work dir goes to /data/scratch/$USER automatically
nextflow run main.nf \
  -profile apocrita,singularity \
  --input    samplesheet.csv \
  --counts   salmon.merged.gene_counts.tsv \
  --contrasts contrasts.csv \
  --organism  human \
  --outdir    /data/projects/myproject/bulkannex_results

# Override the work directory explicitly (all profiles)
nextflow run main.nf \
  -profile apocrita,singularity \
  --work_dir /data/scratch/myuser/nxf_work \
  --input    samplesheet.csv \
  --counts   salmon.merged.gene_counts.tsv \
  --contrasts contrasts.csv \
  --outdir    /data/projects/myproject/bulkannex_results

# Or use Nextflow's native -w flag (equivalent)
nextflow run main.nf -profile apocrita,singularity -w /data/scratch/myuser/nxf_work ...

# Local (Docker)
nextflow run main.nf \
  -profile docker \
  --input    samplesheet.csv \
  --counts   salmon.merged.gene_counts.tsv \
  --contrasts contrasts.csv \
  --outdir    results
```

### Work directory

The `work/` folder holds all intermediate files and can grow very large. Keep it **separate from `--outdir`**:

| Situation | Recommendation |
|-----------|---------------|
| Apocrita HPC | Uses `/data/scratch/$USER/bulkannex_work` by default (via the `apocrita` profile) |
| Other HPC | `--work_dir /scratch/$USER/nxf_work` or `-w /scratch/$USER/nxf_work` |
| Local | Default `./work` is fine; clean with `nextflow clean -f` after successful run |

---

## Launch the dashboard

```bash
bash dashboard/launch_dashboard.sh bulkannex_results 3838
```

Then open `http://localhost:3838` in your browser.

---

## Requirements

- Nextflow >= 23.10
- Singularity, Docker, or conda
- R ≥ 4.3 (for local dashboard only)

---

## Output

See [docs/output.md](docs/output.md) for a full description of output files.


---

## License

MIT
