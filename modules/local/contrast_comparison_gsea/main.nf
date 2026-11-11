process CONTRAST_COMPARISON_GSEA {
    tag "gsea-cross-contrast"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/gsea/contrast_comparison", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    path gsea_files   // all *_GO_all_gsea.csv / *_KEGG_gsea.csv / *_Reactome_gsea.csv

    output:
    path "contrast_comparison_gsea_*.{pdf,png}",          optional: true, emit: plots
    path "contrast_comparison_gsea_common_pathways_*.csv", optional: true, emit: common_pathways
    path 'versions.yml',                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/contrast_comparison_gsea.R \\
        ${args}
    """
}
