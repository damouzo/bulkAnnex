process CONTRAST_COMPARISON {
    tag "cross-contrast"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/dge/contrast_comparison", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    path dge_files   // all *_DESeq2_results.csv collected from all contrasts

    output:
    path "contrast_comparison_*.{pdf,png}", optional: true, emit: plots
    path "contrast_comparison_common_degs.csv",             emit: common_degs
    path 'versions.yml',                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/contrast_comparison.R \\
        ${args}
    """
}
