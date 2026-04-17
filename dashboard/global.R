# global.R — bulkAnnex dashboard global setup
# Sourced once at startup (before ui.R and server.R) by app.R.
# Loads all libraries, resolves the results directory, and defines data-loading helpers.

suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(dplyr)
    library(ggplot2)
    library(plotly)
    library(DT)
    library(enrichplot)
    # ggridges is used by enrichplot::ridgeplot() and loaded as its dependency.
    # The Ridgeplot tab falls back gracefully if it is unavailable.
})

# Null-coalescing operator (rlang-style) — used in modules without requiring rlang on the search path.
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Results directory -------------------------------------------------------
# Priority: BULKANNEX_RESULTS_DIR env var → parent directory of the app → "."
RESULTS_DIR <- local({
    env_val <- Sys.getenv("BULKANNEX_RESULTS_DIR", unset = "")
    if (nchar(env_val) > 0) {
        normalizePath(env_val, mustWork = FALSE)
    } else {
        # Default: parent of the directory containing app.R
        normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
    }
})

message("bulkAnnex dashboard: results dir = ", RESULTS_DIR)

# ---- Data loading functions --------------------------------------------------

load_samplesheet <- function(results_dir) {
    f <- file.path(results_dir, "pipeline_info", "validated_samplesheet.csv")
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    else data.frame()
}

load_vst_counts <- function(results_dir) {
    f <- file.path(results_dir, "normalization", "deseq2_vst_counts.tsv")
    if (file.exists(f)) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    else data.frame()
}

load_qc_metrics <- function(results_dir) {
    f <- file.path(results_dir, "qc", "qc_metrics.csv")
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    else data.frame()
}

load_dge_results <- function(results_dir) {
    dge_dir <- file.path(results_dir, "dge")
    if (!dir.exists(dge_dir)) return(list())
    files <- list.files(dge_dir, pattern = "_DESeq2_results\\.csv$",
                        full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) return(list())
    result <- lapply(files, function(f) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
    names(result) <- sub("_DESeq2_results\\.csv$", "", basename(files))
    result
}

load_gsea_csv <- function(results_dir) {
    gsea_dir <- file.path(results_dir, "gsea")
    if (!dir.exists(gsea_dir)) return(list())
    patterns <- c(
        go       = "_GO_all_gsea\\.csv$",
        kegg     = "_KEGG_gsea\\.csv$",
        reactome = "_Reactome_gsea\\.csv$"
    )
    result <- list()
    for (db in names(patterns)) {
        files <- list.files(gsea_dir, pattern = patterns[[db]],
                            full.names = TRUE, recursive = TRUE)
        if (length(files) > 0) {
            db_list <- lapply(files, function(f) {
                df <- tryCatch(
                    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) {
                        message("Warning: could not parse ", basename(f),
                                " (", e$message, ") — skipping.")
                        NULL
                    }
                )
                if (is.null(df) || nrow(df) == 0) NULL else df
            })
            names(db_list) <- sub(patterns[[db]], "", basename(files))
            result[[db]] <- Filter(Negate(is.null), db_list)
        }
    }
    result
}

# Load native gseaResult / data.table objects saved by gsea_analysis.R.
# Returns list keyed by contrast_id; each value is a named list:
#   list(go_bp, go_mf, go_cc, kegg, reactome)
load_gsea_objects <- function(results_dir) {
    gsea_dir <- file.path(results_dir, "gsea")
    if (!dir.exists(gsea_dir)) return(list())
    rds_files <- list.files(gsea_dir, pattern = "_gsea_dashboard_data\\.rds$",
                            full.names = TRUE, recursive = TRUE)
    if (length(rds_files) == 0) return(list())
    result <- lapply(rds_files, function(f) {
        tryCatch(readRDS(f), error = function(e) {
            message("Warning: could not load ", basename(f), ": ", e$message)
            NULL
        })
    })
    names(result) <- sub("_gsea_dashboard_data\\.rds$", "", basename(rds_files))
    Filter(Negate(is.null), result)
}

# Master loader — returns the full app_data list consumed by all modules.
load_all_data <- function(results_dir = RESULTS_DIR) {
    list(
        samplesheet     = load_samplesheet(results_dir),
        vst             = load_vst_counts(results_dir),
        qc_metrics      = load_qc_metrics(results_dir),
        qc_heatmap_path = file.path(results_dir, "qc", "qc_correlation_heatmap.png"),
        dge             = load_dge_results(results_dir),
        gsea            = load_gsea_csv(results_dir),
        gsea_objects    = load_gsea_objects(results_dir)
    )
}

# Helper: return the path of the pre-rendered plot PNG saved by the pipeline,
# or NULL if it does not exist.
# type: "dotplot" or "ridgeplot"
# Used by mod_gsea to serve fast static images.
get_gsea_png <- function(results_dir, contrast, database, go_ont = "BP",
                         type = "dotplot") {
    gsea_dir <- file.path(results_dir, "gsea", contrast)
    suffix <- if (type == "ridgeplot") "_ridgeplot.png" else "_dotplot.png"
    filename <- switch(database,
        "go"       = paste0(contrast, "_GO_", toupper(go_ont %||% "BP"), suffix),
        "kegg"     = paste0(contrast, "_KEGG", suffix),
        "reactome" = paste0(contrast, "_Reactome", suffix),
        NULL
    )
    if (is.null(filename)) return(NULL)
    f <- file.path(gsea_dir, filename)
    if (file.exists(f)) f else NULL
}

# Helper: return named vector of pathview PNGs for a contrast.
# Names are pathway IDs (e.g., "hsa04110"), values are full file paths.
get_pathview_files <- function(results_dir, contrast) {
    gsea_dir <- file.path(results_dir, "gsea", contrast)
    if (!dir.exists(gsea_dir)) return(character(0))
    files <- list.files(gsea_dir,
                        pattern = paste0("_pathview_.*\\.png$"),
                        full.names = TRUE)
    if (length(files) == 0) return(character(0))
    # Extract pathway IDs: {contrast_id}_pathview_{pathway_id}.png
    pw_ids <- sub(paste0("^", gsub("([][{}()+*?.\\^$|])", "\\\\\\1", contrast),
                         "_pathview_"), "",
                  sub("\\.png$", "", basename(files)))
    names(files) <- pw_ids
    files
}
