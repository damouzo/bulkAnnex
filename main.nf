#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bulkAnnex
    Automated bulkRNA-seq analysis pipeline with interactive Shiny dashboard
    Starting from salmon.merged.gene_counts.tsv produced by nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Homepage : https://github.com/BCI-KRP/bulkAnnex
    Author   : BCI-KRP
    License  : MIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// Print help and exit
if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                          bulkAnnex v1.0                             ║
    ║     Automated bulk RNA-seq downstream analysis + Shiny dashboard     ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf \\
            --input     samplesheet.csv \\
            --counts    salmon.merged.gene_counts.tsv \\
            --contrasts contrasts.csv \\
            --organism  human \\
            --outdir    results

    Mandatory:
        --input         Samplesheet CSV (columns: sample, condition [, batch])
        --counts        Gene counts TSV from nf-core/rnaseq STAR+Salmon
        --contrasts     Contrasts CSV (columns: contrast_id, variable, reference, treatment)
        --organism      Organism: 'human' or 'mouse'

    Optional:
        --outdir            Output directory [results]
        --padj_threshold    Adjusted p-value cutoff [0.05]
        --lfc_threshold     log2FC cutoff for significance [1.0]
        --min_counts        Min count sum per gene for filtering [10]
        --min_samples       Min samples with >= min_counts [2]
        --normalization     Normalization method: vst or rlog [vst]
        --run_gsea          Run GSEA analysis [true]
        --max_memory        Max memory per process [128.GB]
        --max_cpus          Max CPUs per process [16]

    Profiles:
        -profile singularity    Use Singularity (recommended for HPC)
        -profile docker         Use Docker (local)
        -profile conda          Use Conda
        -profile apocrita       QMUL Apocrita HPC (use with singularity)
        -profile slurm          Generic SLURM cluster
        -profile test           Run with demo data
    """.stripIndent()
    exit 0
}

// Validate mandatory parameters
def mandatory = ['input', 'counts', 'contrasts', 'organism']
mandatory.each { param ->
    if (!params[param]) {
        error "Missing mandatory parameter: --${param}. Use --help for usage."
    }
}

if (!['human', 'mouse'].contains(params.organism)) {
    error "Invalid --organism '${params.organism}'. Must be 'human' or 'mouse'."
}

log.info """
╔══════════════════════════════════════════════════════════════════════╗
║                         bulkAnnex v1.0                              ║
╚══════════════════════════════════════════════════════════════════════╝
 input        : ${params.input}
 counts       : ${params.counts}
 contrasts    : ${params.contrasts}
 organism     : ${params.organism}
 outdir       : ${params.outdir}
 padj cutoff  : ${params.padj_threshold}
 lfc cutoff   : ${params.lfc_threshold}
 run_gsea     : ${params.run_gsea}
──────────────────────────────────────────────────────────────────────
""".stripIndent()

// Include main workflow
include { BULKANNEX } from './workflows/bulkannex'

// Entry point
workflow {
    BULKANNEX ()
}

workflow.onComplete {
    log.info (workflow.success ? "\nPipeline completed successfully!\nResults: ${params.outdir}\n" : "\nPipeline failed.\n")
}
