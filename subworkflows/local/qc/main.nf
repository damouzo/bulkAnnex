include { DESEQ2_QC } from '../../../modules/local/deseq2_qc/main'

workflow QC {
    take:
    ch_counts      // path: counts matrix TSV
    ch_samplesheet // path: validated samplesheet CSV

    main:
    ch_versions = Channel.empty()

    DESEQ2_QC(
        ch_counts,
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(DESEQ2_QC.out.versions)

    emit:
    vst_matrix          = DESEQ2_QC.out.vst_matrix
    metrics             = DESEQ2_QC.out.metrics
    library_sizes       = DESEQ2_QC.out.library_sizes
    count_distribution  = DESEQ2_QC.out.count_distribution
    pca                 = DESEQ2_QC.out.pca
    correlation_heatmap = DESEQ2_QC.out.correlation_heatmap
    dispersion          = DESEQ2_QC.out.dispersion
    versions            = ch_versions
}
