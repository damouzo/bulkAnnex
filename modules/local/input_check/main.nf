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
    path 'README.txt',                emit: readme
    path 'versions.yml',              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/validate_inputs.R \\
        --samplesheet ${samplesheet} \\
        --contrasts   ${contrasts} \\
        --counts      ${counts}

    cat > README.txt <<'BULKANNEX_README'
bulkAnnex results
=================

Output layout
-------------
counts/raw/raw_counts.tsv        Rounded integer counts actually passed to DESeq2
                                 (all genes, all samples in the samplesheet).

counts/norm/<norm_group>/        Per-group DESeq2 normalisation:
    deseq2_dds.rds               Fitted DESeqDataSet (reused by DGE).
    deseq2_vst_counts.tsv        VST (blind=FALSE) matrix within the group.
    deseq2_size_factors.csv      Per-sample size factors for the group.

qc/qc_vst_matrix.tsv             Global VST (blind=TRUE) over ALL samples.
                                 Use for cross-group exploration (PCA,
                                 correlation, gene explorer) and QC.

qc/                              QC plots and metrics.
dge/<contrast_id>/               Per-contrast DESeq2 results and plots.
gsea/<contrast_id>/              Per-contrast GSEA results and plots.
pipeline_info/                   Validated inputs, reports, software versions.

Notes
-----
- counts/norm/<norm_group>/deseq2_vst_counts.tsv is design-aware (blind=FALSE)
  and computed within one norm_group; its values are NOT comparable across
  groups, nor with qc/qc_vst_matrix.tsv.
- qc/qc_vst_matrix.tsv is blind (unbiased) and covers all samples; it is the
  recommended matrix for visualisation.
BULKANNEX_README
    """
}
