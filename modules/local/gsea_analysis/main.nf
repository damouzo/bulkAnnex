process GSEA_ANALYSIS {
    tag "${meta.contrast_id}"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/gsea/${meta.contrast_id}", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    tuple val(meta), path(dge_results)

    output:
    tuple val(meta), path("${meta.contrast_id}_GO_all_gsea.csv"),      optional: true, emit: go_results
    tuple val(meta), path("${meta.contrast_id}_KEGG_gsea.csv"),         optional: true, emit: kegg_results
    tuple val(meta), path("${meta.contrast_id}_Hallmarks_gsea.csv"),    optional: true, emit: hallmarks_results
    tuple val(meta), path("${meta.contrast_id}_Reactome_gsea.csv"),     optional: true, emit: reactome_results
    tuple val(meta), path("${meta.contrast_id}_gsea_dashboard_data.rds"), optional: true, emit: dashboard_rds
    path "${meta.contrast_id}_GO_*_dotplot.{pdf,png}",                  optional: true, emit: go_plots
    path "${meta.contrast_id}_KEGG_dotplot.{pdf,png}",                  optional: true, emit: kegg_plots
    path "${meta.contrast_id}_Hallmarks_dotplot.{pdf,png}",             optional: true, emit: hallmarks_plots
    path "${meta.contrast_id}_Reactome_dotplot.{pdf,png}",              optional: true, emit: reactome_plots
    path 'versions.yml',                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/gsea_analysis.R \\
        --dge_results ${dge_results} \\
        --contrast_id ${meta.contrast_id} \\
        --organism    ${params.organism} \\
        ${args}
    """
}
