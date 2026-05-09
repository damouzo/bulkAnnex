process DESEQ2_NORMALIZATION {
    tag "${meta.norm_group}"
    label 'process_medium'

    container params.container

    publishDir [
        path: { meta.norm_group == "all"
            ? "${params.outdir}/normalization"
            : "${params.outdir}/normalization/${meta.norm_group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    ]

    input:
    tuple val(meta), path(samplesheet), path(counts)

    output:
    tuple val(meta), path('deseq2_dds.rds'),          emit: dds
    tuple val(meta), path('deseq2_vst_counts.tsv'),   emit: vst_counts
    tuple val(meta), path('deseq2_size_factors.csv'), emit: size_factors
    path 'versions.yml',                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/normalization.R \\
        --counts      ${counts} \\
        --samplesheet ${samplesheet} \\
        --norm_group  ${meta.norm_group} \\
        --prefix      deseq2 \\
        ${args}
    """
}
