include { DESEQ2_DGE    } from '../../../modules/local/deseq2_dge/main'
include { GSEA_ANALYSIS } from '../../../modules/local/gsea_analysis/main'

workflow DGE_GSEA {
    take:
    ch_contrasts   // channel: [ meta(contrast_id, variable, reference, treatment), dds, vst_counts ]
    ch_run_gsea    // val: boolean — whether to run GSEA

    main:
    ch_versions = Channel.empty()

    // ---- DGE per contrast --------------------------------------------------
    DESEQ2_DGE( ch_contrasts )
    ch_versions = ch_versions.mix(DESEQ2_DGE.out.versions.first())

    // ---- GSEA per contrast (optional) --------------------------------------
    ch_gsea_go        = Channel.empty()
    ch_gsea_kegg      = Channel.empty()
    ch_gsea_hallmarks = Channel.empty()
    ch_gsea_reactome  = Channel.empty()

    if (ch_run_gsea) {
        GSEA_ANALYSIS( DESEQ2_DGE.out.results )
        ch_versions       = ch_versions.mix(GSEA_ANALYSIS.out.versions.first())
        ch_gsea_go        = GSEA_ANALYSIS.out.go_results
        ch_gsea_kegg      = GSEA_ANALYSIS.out.kegg_results
        ch_gsea_hallmarks = GSEA_ANALYSIS.out.hallmarks_results
        ch_gsea_reactome  = GSEA_ANALYSIS.out.reactome_results
    }

    emit:
    dge_results       = DESEQ2_DGE.out.results
    gsea_go           = ch_gsea_go
    gsea_kegg         = ch_gsea_kegg
    gsea_hallmarks    = ch_gsea_hallmarks
    gsea_reactome     = ch_gsea_reactome
    versions          = ch_versions
}
