process DESEQ2_DGE {
    tag "${meta.contrast_id}"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/dge/${meta.contrast_id}", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    tuple val(meta), path(dds), path(vst_counts)

    output:
    tuple val(meta), path("${meta.contrast_id}_DESeq2_results.csv"), emit: results
    path "${meta.contrast_id}_volcano.{pdf,png}",                    emit: volcano
    path "${meta.contrast_id}_MA.{pdf,png}",                         emit: ma_plot
    path "${meta.contrast_id}_top_genes_heatmap.{pdf,png}",          optional: true, emit: heatmap
    path 'versions.yml',                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def contrast_spec = "${meta.variable},${meta.reference},${meta.treatment}"
    """
    Rscript ${projectDir}/bin/dge_analysis.R \\
        --dds         ${dds} \\
        --contrast    "${contrast_spec}" \\
        --contrast_id ${meta.contrast_id} \\
        --vst         ${vst_counts} \\
        ${args}
    """
}
