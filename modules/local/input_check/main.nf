process INPUT_CHECK {
    tag "input validation"
    label 'process_single'

    container params.container

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename == 'versions.yml') null else filename
    }

    input:
    path samplesheet
    path contrasts
    path counts

    output:
    path 'validated_samplesheet.csv', emit: samplesheet
    path 'validated_contrasts.csv',   emit: contrasts
    path 'versions.yml',              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/validate_inputs.R \\
        --samplesheet ${samplesheet} \\
        --contrasts   ${contrasts} \\
        --counts      ${counts}
    """
}
