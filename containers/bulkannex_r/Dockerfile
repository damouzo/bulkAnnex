FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="bulkannex_r"
LABEL org.opencontainers.image.description="R/Bioconductor environment for bulkAnnex bulk RNA-seq downstream analysis pipeline"
LABEL org.opencontainers.image.version="1.0.1"
LABEL org.opencontainers.image.licenses="MIT"

# ---- System dependencies -----------------------------------------------------
# libglpk-dev + libgmp3-dev: required by igraph (clusterProfiler dep).
# All other system libs (libcurl, libssl, libxml2, libcairo, libpango)
# are already present in bioconductor_docker:RELEASE_3_20.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libglpk-dev \
        libgmp3-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- Bioconductor packages ---------------------------------------------------
# Installed via BiocManager — guarantees the FULL (non-pruned) DESeq2 including
# GenomeInfoDbData, which bioconda's bioconductor-deseq2_pruned actively excludes.
# bioconductor_docker:RELEASE_3_20 ships BiocManager pre-configured for release 3.20.
RUN Rscript -e " \
    BiocManager::install( \
        c( \
            'GenomeInfoDbData', \
            'DESeq2', \
            'BiocParallel', \
            'clusterProfiler', \
            'ReactomePA', \
            'fgsea', \
            'org.Hs.eg.db', \
            'org.Mm.eg.db', \
            'AnnotationDbi', \
            'pathview', \
            'edgeR' \
        ), \
        ask    = FALSE, \
        update = FALSE \
    ) \
"

# ---- CRAN packages -----------------------------------------------------------
RUN Rscript -e " \
    install.packages( \
        c( \
            'optparse', \
            'ashr', \
            'ggplot2', \
            'ggrepel', \
            'ggridges', \
            'pheatmap', \
            'RColorBrewer', \
            'viridis', \
            'patchwork', \
            'scales', \
            'forcats', \
            'dplyr', \
            'tidyr', \
            'readr', \
            'tibble', \
            'stringr', \
            'jsonlite', \
            'rlang', \
            'shiny', \
            'bslib', \
            'DT', \
            'plotly', \
            'htmlwidgets', \
            'shinyWidgets', \
            'shinycssloaders' \
        ), \
        repos = 'https://cloud.r-project.org' \
    ) \
"

# ---- Sanity check ------------------------------------------------------------
RUN Rscript -e " \
    pkgs <- c('DESeq2', 'GenomeInfoDbData', 'clusterProfiler', 'ReactomePA', 'fgsea', \
              'pathview', 'org.Hs.eg.db', 'org.Mm.eg.db', 'ashr', \
              'ggridges', 'shiny'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Missing packages: ', paste(missing, collapse=', ')); \
    message('All required packages present.') \
"
