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
| `--input` | `samplesheet.csv` — columns: `sample`, `condition`, `batch` (optional) |
| `--counts` | `salmon.merged.gene_counts.tsv` — gene_id, gene_name, one column per sample |
| `--contrasts` | `contrasts.csv` — columns: `contrast_id`, `variable`, `reference`, `treatment` |
| `--organism` | `human` (default) or `mouse` |
| `--outdir` | Output directory (default: `results`) |

---

## Pipeline steps

1. **Input validation** — checks samplesheet, contrasts, and counts consistency
2. **QC** — library sizes, count distributions, PCA, correlation heatmap, dispersion
3. **Normalization** — DESeq2 VST, fitted DDS saved for DGE reuse
4. **DGE** — DESeq2 results per contrast, volcano, MA, top-gene heatmap
5. **GSEA** — GO (BP/MF/CC), KEGG, MSigDB Hallmarks, Reactome via clusterProfiler/fgsea
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
