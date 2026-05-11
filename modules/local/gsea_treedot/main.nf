process GSEA_TREEDOT {
    tag "cross-contrast"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/gsea/treedot", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    path dge_files   // all *_DESeq2_results.csv collected from all contrasts

    output:
    path "gsea_treedot_*.{pdf,png}", optional: true, emit: plots
    path "gsea_treedot_*.rds",       optional: true, emit: rds
    path 'versions.yml',              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/gsea_treedot.R \\
        --organism ${params.organism} \\
        ${args}
    """
}
