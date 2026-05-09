# Usage

## Samplesheet format

```csv
sample,condition,batch,norm_group
NB4_Exp4_Scramble,NB4_Scramble,Exp4,NB4
NB4_Exp4_DDX41sh1,NB4_DDX41_sh1,Exp4,NB4
MSCline_Exp8_Scramble,MSCline_Scramble,Exp8,MSCline
DDX41Patient_MNC_P1,MNC_DDX41_Patient,,MNC
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample name. Must exactly match a column header in the counts matrix. No spaces. |
| `condition` | Yes | Experimental condition. Must be a valid R identifier (use underscores, not spaces or hyphens). |
| `batch` | No | Batch covariate. If present, DESeq2 design becomes `~ batch + condition`. All rows must have a value if the column exists. |
| `norm_group` | No | Normalisation group. Samples sharing the same value are fitted in a single DESeq2 model. Omit (or set to `all`) to normalise all samples together. See [Group-aware normalisation](#group-aware-normalisation) below. |

## Contrasts format

```csv
contrast_id,variable,reference,treatment,norm_group
NB4_DDX41_sh1_vs_Scramble,condition,NB4_Scramble,NB4_DDX41_sh1,NB4
MSCline_DDX41_sh1_vs_Scramble,condition,MSCline_Scramble,MSCline_DDX41_sh1,MSCline
MNC_DDX41_Patient_vs_Healthy,condition,MNC_Healthy,MNC_DDX41_Patient,MNC
```

| Column | Required | Description |
|--------|----------|-------------|
| `contrast_id` | Yes | Unique identifier for the contrast. Used as output directory name and file prefix. No spaces. |
| `variable` | Yes | Column in the samplesheet to test (typically `condition`). |
| `reference` | Yes | Reference level (denominator, e.g. control). Must be a value present in the `variable` column **within the norm_group**. |
| `treatment` | Yes | Treatment level (numerator). Must be a value present in the `variable` column **within the norm_group**. |
| `norm_group` | No | Must match a `norm_group` value in the samplesheet. The contrast is computed against the DESeq2 model fitted for that group. Omit to use the global model. |

## Group-aware normalisation

When an experiment mixes biologically distinct sample types (e.g. cell lines and primary patient cells), a single global DESeq2 model can distort dispersion estimates. Adding `norm_group` runs one independent DESeq2 normalisation per group:

- Each `norm_group` runs DESeq2 on its own samples only.
- QC (PCA, library sizes, heatmap) always uses **all samples together**.
- `dge/` and `gsea/` output directories remain flat — contrast IDs are globally unique.
- Omitting `norm_group` is fully backward compatible (equivalent to `norm_group = all`).

See `data_demo/samplesheet.csv` and `data_demo/contrasts.csv` for a complete worked example with four independent normalisation groups (NB4, MSCline, MNC, MSC).

## Counts matrix format

Tab-separated, as produced by nf-core/rnaseq:

```
gene_id    gene_name    sample1    sample2    ...
ENSG0001   TP53         1234.5     2345.0     ...
```

Salmon float counts are rounded to integers automatically.

## Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Use Singularity containers |
| `docker` | Use Docker containers |
| `conda` | Use conda environments |
| `apocrita` | Apocrita HPC (SLURM executor) |
| `slurm` | Generic SLURM cluster |
| `local` | Local execution (no job scheduler) |
| `test` | Minimal CI run with demo data |

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--organism` | `human` | `human` or `mouse` |
| `--run_gsea` | `true` | Enable/disable GSEA |
| `--skip_dashboard` | `false` | Skip Shiny dashboard step |
| `--alpha` | `0.05` | DESeq2 FDR threshold |
| `--lfc_threshold` | `0.0` | LFC threshold for DESeq2 results() |
| `--gsea_min_gs` | `15` | Minimum gene set size |
| `--gsea_max_gs` | `500` | Maximum gene set size |
| `--max_memory` | `128.GB` | Resource cap |
| `--max_cpus` | `16` | Resource cap |
| `--max_time` | `240.h` | Resource cap |
| `--work_dir` | `null` | Override work directory (intermediate files). Equivalent to `-w /path`. |

## Work directory

Nextflow writes all intermediate files to a `work/` directory. This can grow large — keep it on a scratch filesystem, not alongside your results.

```bash
# Option 1: use --work_dir param
nextflow run main.nf -profile apocrita,singularity \
  --work_dir /data/scratch/myuser/nxf_work \
  --outdir   /data/projects/myproject/results ...

# Option 2: use Nextflow's native -w flag (identical effect)
nextflow run main.nf -profile apocrita,singularity \
  -w /data/scratch/myuser/nxf_work \
  --outdir /data/projects/myproject/results ...
```

The `apocrita` profile already sets `workDir = "/data/scratch/$USER/bulkannex_work"` as default.  
Clean up after a successful run with: `nextflow clean -f`
