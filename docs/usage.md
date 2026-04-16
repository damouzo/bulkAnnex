# Usage

## Samplesheet format

```csv
sample,condition,batch
NB4_Exp4_Scramble,Scramble,Exp4
NB4_Exp4_DDX41sh1,DDX41_sh1,Exp4
```

- `sample` must exactly match column names in the counts matrix.
- `condition` is required. Must be a valid R variable name (no spaces or hyphens — use underscores).
- `batch` is optional. If present, DESeq2 design becomes `~ batch + condition`.

## Contrasts format

```csv
contrast_id,variable,reference,treatment
DDX41_sh1_vs_Scramble,condition,Scramble,DDX41_sh1
```

- `contrast_id`: used as output file prefix. No spaces.
- `variable`: must match a column in the samplesheet.
- `reference`: the denominator group (e.g., control).
- `treatment`: the numerator group.

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
