process SINGULARITY_PULL {
    tag "pull_container"
    label 'process_medium'   // 4 CPUs, 16 GB — enough for mksquashfs

    // Run natively (no container) on the SLURM compute node.
    // This guarantees mksquashfs has the RAM it needs, unlike the login node.
    container null

    output:
    val(true), emit: ready

    when:
    task.ext.when == null || task.ext.when

    script:
    def sif   = "${projectDir}/containers/bulkannex_r/bulkannex_r_1.0.0.sif"
    def image = "docker://damouzo/bulkannex_r:1.0.0"
    // Use a job-unique scratch dir for mksquashfs temporary files.
    // SLURM_JOB_ID is set in every SLURM job; fall back to PID for local runs.
    def tmpbase = "/gpfs/scratch/\${USER:-tmp}/singularity_tmp"
    """
    SIF="${sif}"
    IMAGE="${image}"
    LOCK="\${SIF}.lock"
    TMPDIR="${tmpbase}_\${SLURM_JOB_ID:-\$\$}"

    # ---- Fast path: SIF already present -------------------------------------
    if [ -f "\${SIF}" ]; then
        echo "SIF already exists — skipping pull."
        exit 0
    fi

    # ---- Advisory lock to serialise concurrent pulls (cross-node via GPFS) --
    WAITED=0
    while [ -f "\${LOCK}" ] && [ \${WAITED} -lt 3600 ]; do
        echo "Another pull in progress (lock: \${LOCK}). Waiting 30 s..."
        sleep 30
        WAITED=\$((WAITED + 30))
    done

    # Re-check: a concurrent job may have finished while we waited
    if [ -f "\${SIF}" ]; then
        echo "SIF appeared while waiting — skipping pull."
        exit 0
    fi

    # ---- Pull ---------------------------------------------------------------
    mkdir -p "\${TMPDIR}"
    export SINGULARITY_TMPDIR="\${TMPDIR}"
    export APPTAINER_TMPDIR="\${TMPDIR}"

    touch "\${LOCK}"
    trap 'rm -f "\${LOCK}"; rm -rf "\${TMPDIR}"' EXIT

    echo "Pulling \${IMAGE} → \${SIF}"
    echo "Node: \$(hostname)"
    echo "Tmp : \${TMPDIR}"
    singularity pull --force "\${SIF}" "\${IMAGE}"

    echo "Pull complete: \${SIF}"
    """
}
