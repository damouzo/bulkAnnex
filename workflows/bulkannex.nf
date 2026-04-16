/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BULKANNEX MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downstream bulk RNA-seq analysis: QC → normalization → DGE → GSEA → Dashboard
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK           } from '../modules/local/input_check/main'
include { DESEQ2_NORMALIZATION  } from '../modules/local/deseq2_normalization/main'
include { SHINY_DASHBOARD       } from '../modules/local/shiny_dashboard/main'
include { SINGULARITY_PULL      } from '../modules/local/singularity_pull/main'
include { QC                    } from '../subworkflows/local/qc/main'
include { DGE_GSEA              } from '../subworkflows/local/dge_gsea/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: BULKANNEX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKANNEX {
    ch_versions = Channel.empty()

    // ---- Stage inputs -------------------------------------------------------
    ch_samplesheet = Channel.fromPath(params.input,     checkIfExists: true)
    ch_contrasts   = Channel.fromPath(params.contrasts, checkIfExists: true)
    ch_counts      = Channel.fromPath(params.counts,    checkIfExists: true)

    // ---- 0. Ensure the Singularity SIF exists on a compute node -------------
    // SINGULARITY_PULL runs natively (no container) on SLURM so that mksquashfs
    // has real RAM. It skips immediately if the SIF already exists.
    // Evaluated at parse time: if SIF is present from the start, skip entirely.
    def sif_path = "${projectDir}/containers/bulkannex_r/bulkannex_r_1.0.0.sif"
    ch_container_ready = Channel.value(true)  // default: SIF exists, proceed immediately

    if (workflow.containerEngine == 'singularity' && !new File(sif_path).exists()) {
        SINGULARITY_PULL()
        ch_container_ready = SINGULARITY_PULL.out.ready
    }

    // Inject the ready signal into the first process (INPUT_CHECK) so that
    // all downstream processes are implicitly gated behind it.
    // combine() creates the dependency; map{} strips the ready value back to a Path.
    ch_samplesheet = ch_samplesheet
        .combine(ch_container_ready)
        .map { ss, _ready -> ss }

    // ---- 1. Validate inputs -------------------------------------------------
    INPUT_CHECK(
        ch_samplesheet,
        ch_contrasts,
        ch_counts
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // ---- 2. QC (raw counts + blind VST) -------------------------------------
    QC(
        ch_counts,
        INPUT_CHECK.out.samplesheet
    )
    ch_versions = ch_versions.mix(QC.out.versions)

    // ---- 3. Normalization (fit full DESeq2 model, save DDS) -----------------
    DESEQ2_NORMALIZATION(
        ch_counts,
        INPUT_CHECK.out.samplesheet
    )
    ch_versions = ch_versions.mix(DESEQ2_NORMALIZATION.out.versions)

    // ---- 4. Build per-contrast channel with meta map -----------------------
    //  Read contrasts CSV → one row per contrast → combine with DDS + VST
    ch_contrast_rows = INPUT_CHECK.out.contrasts
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                contrast_id : row.contrast_id,
                variable    : row.variable,
                reference   : row.reference,
                treatment   : row.treatment
            ]
            return meta
        }

    // Combine each contrast meta with the shared DDS and VST outputs
    ch_dge_input = ch_contrast_rows
        .combine(DESEQ2_NORMALIZATION.out.dds)
        .combine(DESEQ2_NORMALIZATION.out.vst_counts)
        .map { meta, dds, vst -> [ meta, dds, vst ] }

    // ---- 5. DGE + GSEA per contrast ----------------------------------------
    DGE_GSEA(
        ch_dge_input,
        params.run_gsea
    )
    ch_versions = ch_versions.mix(DGE_GSEA.out.versions)

    // ---- 6. Shiny dashboard -------------------------------------------------
    // The dashboard reads results at runtime from outdir/ via BULKANNEX_RESULTS_DIR.
    // We only need a sync signal (collected versions) so the process runs after
    // all upstream processes have published their outputs.
    if (!params.skip_dashboard) {
        ch_sync = ch_versions.collect()
        SHINY_DASHBOARD(ch_sync)
        ch_versions = ch_versions.mix(SHINY_DASHBOARD.out.versions)
    }

    // ---- 7. Pipeline info ---------------------------------------------------
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}
