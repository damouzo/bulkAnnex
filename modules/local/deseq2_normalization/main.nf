process DESEQ2_NORMALIZATION {
    tag "normalization"
    label 'process_medium'

    container params.container

    publishDir "${params.outdir}/normalization", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    path counts
    path samplesheet

    output:
    path 'deseq2_dds.rds',          emit: dds
    path 'deseq2_vst_counts.tsv',   emit: vst_counts
    path 'deseq2_size_factors.csv', emit: size_factors
    path 'versions.yml',            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/normalization.R \\
        --counts      ${counts} \\
        --samplesheet ${samplesheet} \\
        --prefix      deseq2 \\
        ${args}
    """
}
