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
    # ggridges is used internally by enrichplot::ridgeplot(); not loaded here
    # so a missing package does not crash the app — ridgeplot tab shows an error message instead.
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
        go        = "_GO_all_gsea\\.csv$",
        kegg      = "_KEGG_gsea\\.csv$",
        hallmarks = "_Hallmarks_gsea\\.csv$",
        reactome  = "_Reactome_gsea\\.csv$"
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
#   list(go_bp, go_mf, go_cc, kegg, hallmarks, reactome)
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
