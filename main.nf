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
    if (workflow.success) {
        // ---- Always refresh dashboard scripts from projectDir ----------------
        // The SHINY_DASHBOARD process is cached with -resume; re-publishing from
        // the old work dir would overwrite any script updates. Copying here in
        // onComplete (which always runs) ensures the latest scripts reach outdir.
        def dashSrc = new File("${projectDir}/dashboard")
        def dashDst = new File("${params.outdir}/dashboard")
        dashDst.mkdirs()
        new File(dashDst, 'modules').mkdirs()
        new File(dashDst, 'www').mkdirs()

        // Top-level dashboard files
        ['app.R', 'global.R', 'ui.R', 'server.R',
         'launch_dashboard.sh', 'launch_dashboard_hpc.sh'].each { fname ->
            def src = new File(dashSrc, fname)
            def dst = new File(dashDst, fname)
            if (src.exists()) { dst.bytes = src.bytes }
        }
        // modules/ and www/ subdirectories
        ['modules', 'www'].each { subdir ->
            def srcSubDir = new File(dashSrc, subdir)
            if (srcSubDir.exists()) {
                srcSubDir.eachFile { f ->
                    new File(dashDst, "${subdir}/${f.name}").bytes = f.bytes
                }
            }
        }

        // Resolve SIF path so the HPC command is copy-paste ready (no manual editing).
        def sif_abs  = new File("${projectDir}/containers/bulkannex_r/bulkannex_r_1.0.0.sif").canonicalPath
        def sif_line = new File(sif_abs).exists()
            ? "BULKANNEX_SIF=${sif_abs} \\\n     bash"
            : "bash"
        log.info """
──────────────────────────────────────────────────────────────────────
 Pipeline completed successfully!

 Results : ${params.outdir}

 Launch the interactive dashboard
 ─────────────────────────────────────────────────────────────────────
 Local (R must be available):
   bash ${params.outdir}/dashboard/launch_dashboard.sh

 HPC (Apocrita / SLURM) — copy-paste this command:
   ${sif_line} ${params.outdir}/dashboard/launch_dashboard_hpc.sh ${params.outdir}

   Once the SLURM job starts, the exact ssh tunnel command will be
   printed. Run it in a new terminal on your local machine.
──────────────────────────────────────────────────────────────────────
""".stripIndent()
    } else {
        log.info "\nPipeline failed. Check the error above and rerun with -resume.\n"
    }
}
