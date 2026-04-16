process DESEQ2_QC {
    tag "QC"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    path counts
    path samplesheet

    output:
    path 'qc_vst_matrix.tsv',              emit: vst_matrix
    path 'qc_metrics.csv',                 emit: metrics
    path 'qc_library_sizes.{pdf,png}',     emit: library_sizes
    path 'qc_count_distribution.{pdf,png}',emit: count_distribution
    path 'qc_pca.{pdf,png}',              emit: pca
    path 'qc_correlation_heatmap.{pdf,png}', emit: correlation_heatmap
    path 'qc_dispersion.{pdf,png}',        emit: dispersion
    path 'versions.yml',                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/qc_analysis.R \\
        --counts      ${counts} \\
        --samplesheet ${samplesheet} \\
        --prefix      qc \\
        ${args}
    """
}
