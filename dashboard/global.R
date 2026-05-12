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
    library(ggtree)
    library(cowplot)
    library(clusterProfiler)
    library(RColorBrewer)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
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
    # Always prefer the global blind VST from QC (all samples, blind=TRUE).
    # This is the standard DESeq2 approach for visualisation (PCA, Gene Explorer,
    # Count Distribution) — it is not confounded by norm_group boundaries and is
    # 100% valid for cross-group comparisons in exploratory plots.
    #
    # The per-norm_group VSTs in normalization/<group>/deseq2_vst_counts.tsv are
    # only needed for DGE; when multiple groups run they OVERWRITE each other
    # under normalization/ (flat publish), so we NEVER use that file for the
    # dashboard.
    qc_vst <- file.path(results_dir, "qc", "qc_vst_matrix.tsv")
    if (file.exists(qc_vst)) {
        return(read.delim(qc_vst, check.names = FALSE, stringsAsFactors = FALSE))
    }
    # Legacy single-group runs that pre-date the QC blind-VST output.
    flat <- file.path(results_dir, "normalization", "deseq2_vst_counts.tsv")
    if (file.exists(flat)) {
        return(read.delim(flat, check.names = FALSE, stringsAsFactors = FALSE))
    }
    data.frame()
}

load_qc_metrics <- function(results_dir) {
    f <- file.path(results_dir, "qc", "qc_metrics.csv")
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    else data.frame()
}

# Returns validated_contrasts.csv with norm_group column (always present after
# pipeline run; falls back to empty data.frame for pre-norm_group results).
load_contrasts_meta <- function(results_dir) {
    f <- file.path(results_dir, "pipeline_info", "validated_contrasts.csv")
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
        qc_basedir      = file.path(results_dir, "qc"),
        qc_heatmap_path = file.path(results_dir, "qc", "qc_correlation_heatmap.png"),
        dge             = load_dge_results(results_dir),
        gsea            = load_gsea_csv(results_dir),
        gsea_objects    = load_gsea_objects(results_dir),
        contrasts_meta  = load_contrasts_meta(results_dir),
        treedot_rds     = load_treedot_rds(results_dir)
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

# Helper: pre-rendered TreeDot PNG saved by the GSEA_TREEDOT pipeline process.
# tag: e.g. "GO_BP", "KEGG", "Reactome". If NULL, tries GO_BP then legacy file.
get_treedot_png <- function(results_dir, tag = NULL) {
    td_dir <- file.path(results_dir, "gsea", "treedot")
    if (!is.null(tag)) {
        tag_upper <- toupper(tag)
        candidates <- if (tag_upper == "KEGG") {
            c("KEGG", "Kegg")
        } else if (tag_upper == "REACTOME") {
            c("Reactome", "reactome")
        } else if (grepl("^GO_", tag_upper)) {
            c(tag_upper)
        } else {
            c(tag)
        }
        for (cand in unique(candidates)) {
            f <- file.path(td_dir, paste0("gsea_treedot_", cand, ".png"))
            if (file.exists(f)) return(f)
        }
        return(NULL)
    }
    # No tag: try GO_BP first, then legacy generic file
    f <- file.path(td_dir, "gsea_treedot_GO_BP.png")
    if (file.exists(f)) return(f)
    f <- file.path(td_dir, "gsea_treedot.png")
    if (file.exists(f)) return(f)
    NULL
}

# Helper: load all compareClusterResult RDS files saved by GSEA_TREEDOT.
# Returns a named list, e.g. list(GO_BP = <cmp>, KEGG = <cmp>, Reactome = <cmp>).
# Backward-compatible: if only the legacy gsea_treedot.rds exists, returns
# list(GO_BP = <cmp>) so the dashboard still finds it under the GO_BP key.
load_treedot_rds <- function(results_dir) {
    td_dir   <- file.path(results_dir, "gsea", "treedot")
    tagged   <- Sys.glob(file.path(td_dir, "gsea_treedot_*.rds"))
    out_list <- list()
    for (f in sort(tagged)) {
        raw_tag <- sub("^gsea_treedot_(.+)\\.rds$", "\\1", basename(f))
        tag_upper <- toupper(raw_tag)
        tag <- if (tag_upper == "KEGG") {
            "KEGG"
        } else if (tag_upper == "REACTOME") {
            "Reactome"
        } else if (grepl("^GO_", tag_upper)) {
            tag_upper
        } else {
            raw_tag
        }
        cmp <- tryCatch(readRDS(f), error = function(e) {
            message("Warning: could not load ", basename(f), ": ", e$message)
            NULL
        })
        if (!is.null(cmp)) out_list[[tag]] <- cmp
    }
    if (length(out_list) > 0) return(out_list)
    # Backward compat: legacy single-file results
    legacy <- file.path(td_dir, "gsea_treedot.rds")
    if (file.exists(legacy)) {
        cmp <- tryCatch(readRDS(legacy), error = function(e) NULL)
        if (!is.null(cmp)) return(list(GO_BP = cmp))
    }
    NULL
}

# plot_treedot — TreeDot plot from a compareClusterResult object.
# Used both by bin/gsea_treedot.R (pipeline) and mod_gsea.R (dashboard interactive).
# cmp      : compareClusterResult with pairwise_termsim() already applied
# keytype  : gene ID type in core_enrichment ("ENSEMBL" for GO, "ENTREZID" for KEGG)
plot_treedot <- function(cmp, top_paths = 5, clust_num = 3,
                         ora_type = "GO", ora_ont = "BP",
                         ora_min_gs = 10, ora_max_gs = 500, ora_padj = 1,
                         org_db_str = "org.Hs.eg.db", kegg_org = "hsa",
                         keytype = "ENSEMBL") {

    fort <- fortify(cmp, showCategory = top_paths, includeAll = TRUE, split = NULL)
    fort$Cluster <- sub("\n.*", "", fort$Cluster)
    fort$geneID  <- fort$core_enrichment
    if (nrow(fort) == 0) stop("No pathways after fortify.")

    id_unique <- unique(fort$Description)
    cl_unique <- unique(fort$Cluster)
    mat <- matrix(0, nrow = length(id_unique), ncol = length(cl_unique),
                  dimnames = list(id_unique, cl_unique))
    for (i in seq_len(nrow(fort))) mat[fort$Description[i], fort$Cluster[i]] <- 1
    id_mat <- as.data.frame(mat)

    fill_termsim <- function(x, keep) {
        ts <- x@termsim[keep, keep, drop = FALSE]
        ts[is.na(ts)] <- 0
        ts2 <- ts + t(ts)
        diag(ts2) <- 1
        ts2
    }
    ts2 <- tryCatch(
        fill_termsim(cmp, rownames(id_mat)),
        error = function(e) {
            diag_m <- diag(nrow(id_mat))
            rownames(diag_m) <- colnames(diag_m) <- rownames(id_mat)
            diag_m
        }
    )
    hc        <- stats::hclust(stats::as.dist(1 - ts2), method = "ward.D")
    clust_num <- min(clust_num, nrow(id_mat))
    clus      <- stats::cutree(hc, clust_num)

    keywords <- character(clust_num)
    for (i in seq_len(clust_num)) {
        paths_i    <- names(clus[clus == i])
        cluster_df <- fort[fort$Description %in% paths_i, ]
        genes_i    <- unique(unlist(stringr::str_split(cluster_df$geneID, "/")))
        genes_i    <- genes_i[nzchar(genes_i) & !is.na(genes_i)]

        # KEGG and Reactome ORA require ENTREZ IDs; convert if GSEA used ENSEMBL keys
        genes_for_ora <- genes_i
        if (ora_type %in% c("KEGG", "Reactome") && keytype != "ENTREZID") {
            genes_for_ora <- tryCatch(
                na.omit(as.character(AnnotationDbi::mapIds(
                    get(org_db_str), keys = genes_i,
                    column = "ENTREZID", keytype = keytype,
                    multiVals = "first"))),
                error = function(e) character(0)
            )
        }

        reactome_org_ora <- if (grepl("Mm", org_db_str)) "mouse" else "human"
        ora_res <- tryCatch({
            if (ora_type == "GO") {
                clusterProfiler::enrichGO(
                    gene = genes_for_ora, OrgDb = org_db_str, keyType = keytype,
                    ont = ora_ont, pvalueCutoff = ora_padj, pAdjustMethod = "BH",
                    qvalueCutoff = 1, minGSSize = ora_min_gs, maxGSSize = ora_max_gs)
            } else if (ora_type == "KEGG") {
                clusterProfiler::enrichKEGG(
                    gene = genes_for_ora, organism = kegg_org,
                    pvalueCutoff = ora_padj, pAdjustMethod = "BH",
                    qvalueCutoff = 1, minGSSize = ora_min_gs, maxGSSize = ora_max_gs)
            } else {
                ReactomePA::enrichPathway(
                    gene = genes_for_ora, organism = reactome_org_ora,
                    pvalueCutoff = ora_padj, pAdjustMethod = "BH",
                    qvalueCutoff = 1, minGSSize = ora_min_gs, maxGSSize = ora_max_gs)
            }
        }, error = function(e) NULL)

        keywords[i] <- if (!is.null(ora_res) && nrow(ora_res@result) > 0 &&
                            !is.na(ora_res@result$Description[1]))
            ora_res@result$Description[1]
        else
            paste0("Cluster ", i)
    }

    kw_maxlen <- max(nchar(keywords), 17)

    g_split    <- split(names(clus), clus)
    p_tree     <- ggtree::ggtree(hc, size = 1.2)
    clades     <- sapply(g_split, function(n) ggtree::MRCA(p_tree, n))
    p_tree     <- ggtree::groupClade(p_tree, clades, group_name = "SubTree_ORA") +
        ggplot2::aes(color = SubTree_ORA)
    pal        <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(6, "Set1"))

    ggtree_full  <- p_tree +
        ggplot2::scale_color_manual(name = "SubTree ORA",
                                    values = pal, breaks = seq_len(clust_num), labels = keywords) +
        ggplot2::theme(legend.position = "right", legend.justification = c(0, 1.5))
    ggtree_noleg <- ggtree_full + ggplot2::theme(legend.position = "none")

    fort$log_p.adjust <- -log10(fort$p.adjust)
    fort$Description  <- factor(fort$Description, levels = hc$labels[hc$order])
    cl_levels <- if (is.factor(cmp@compareClusterResult$Cluster))
        sub("\n.*", "", levels(cmp@compareClusterResult$Cluster))
    else
        sub("\n.*", "", unique(as.character(cmp@compareClusterResult$Cluster)))
    fort$Cluster <- factor(fort$Cluster, levels = cl_levels)

    dotplot_full <- fort %>%
        ggplot2::ggplot(ggplot2::aes(x = Cluster, y = Description,
                                     color = NES, size = log_p.adjust)) +
        ggplot2::geom_point() +
        ggplot2::scale_y_discrete(position = "right") +
        ggplot2::scale_color_gradient2(low = "blue4", mid = "white", high = "red") +
        cowplot::theme_cowplot(font_size = 15) +
        ggplot2::theme(axis.line   = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
                       axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::xlab("Contrast") + ggplot2::ylab("") +
        ggplot2::guides(size  = ggplot2::guide_legend(title = "-log10(p.adj)"),
                        color = ggplot2::guide_colorbar(title = "NES")) +
        ggplot2::theme(axis.ticks    = ggplot2::element_blank(),
                       legend.position  = "right", legend.justification = c(0, 0),
                       legend.text      = ggplot2::element_text(size = 13),
                       legend.title     = ggplot2::element_text(size = 14))
    dotplot_noleg <- dotplot_full + ggplot2::theme(legend.position = "none")

    leg_tree <- cowplot::get_legend(ggtree_full)
    leg_dot  <- cowplot::get_legend(dotplot_full)
    row_main <- cowplot::plot_grid(ggtree_noleg, NULL, dotplot_noleg,
                                   nrow = 1, rel_widths = c(0.3, -0.05, 2), align = "h")
    leg_col  <- cowplot::plot_grid(leg_dot, NULL, leg_tree,
                                   ncol = 1, rel_heights = c(1, -0.5, 1))
    cowplot::plot_grid(row_main, leg_col, NULL,
                       nrow = 1, rel_widths = c(1, 0.1, kw_maxlen * 0.006))
}

# Helper: extract database info from a compareClusterResult object.
# Returns list(database, go_ont, keytype) — used by the TreeDot server
# to decide whether to reuse the stored RDS or recompute from DGE data.
get_rds_db_info <- function(cmp) {
    default <- list(database = "GO", go_ont = "BP", keytype = "ENSEMBL")
    if (is.null(cmp)) return(default)
    fun_name <- tryCatch(as.character(cmp@.call$fun), error = function(e) "gseGO")
    if (grepl("Pathway|Reactome", fun_name, ignore.case = TRUE))
        return(list(database = "Reactome", go_ont = NA_character_, keytype = "ENTREZID"))
    if (grepl("KEGG", fun_name, ignore.case = TRUE))
        return(list(database = "KEGG", go_ont = NA_character_, keytype = "ENTREZID"))
    go_ont  <- tryCatch(as.character(cmp@.call$ont  %||% "BP"),     error = function(e) "BP")
    keytype <- tryCatch(as.character(cmp@.call$keyType %||% "ENSEMBL"), error = function(e) "ENSEMBL")
    list(database = "GO", go_ont = go_ont, keytype = keytype)
}

# Helper: rebuild a compareClusterResult from DGE data frames for interactive recompute.
# dge_list   : named list of DGE data.frames (from data()$dge) with columns gene_id,
#              log2FoldChange, pvalue
# org_db_str : "org.Hs.eg.db" or "org.Mm.eg.db"
# kegg_org   : "hsa" or "mmu"
# database   : "GO" or "KEGG"
# go_ont     : "BP", "MF", or "CC" (ignored for KEGG)
build_treedot_cmp <- function(dge_list, org_db_str, kegg_org,
                               database = "GO", go_ont = "BP",
                               min_gs = 15L, max_gs = 500L) {
    org_db <- get(org_db_str)
    ranked_ens <- list()
    ranked_ez  <- list()

    for (nm in names(dge_list)) {
        df <- dge_list[[nm]]
        req_cols <- c("gene_id", "log2FoldChange", "pvalue")
        if (!all(req_cols %in% colnames(df))) next
        df <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), ]
        if (nrow(df) == 0) next
        df$rank_metric <- sign(df$log2FoldChange) * (-log10(df$pvalue + .Machine$double.eps))

        rl <- setNames(df$rank_metric, df$gene_id)
        rl <- sort(rl, decreasing = TRUE)
        ranked_ens[[nm]] <- rl[!duplicated(names(rl))]

        ez <- tryCatch(
            AnnotationDbi::mapIds(org_db, keys = df$gene_id, column = "ENTREZID",
                                  keytype = "ENSEMBL", multiVals = "first"),
            error = function(e) NULL
        )
        if (!is.null(ez)) {
            df$entrez_id <- ez[df$gene_id]
            df_ez <- df[!is.na(df$entrez_id), ]
            if (nrow(df_ez) > 0) {
                rl_ez <- setNames(df_ez$rank_metric, df_ez$entrez_id)
                rl_ez <- sort(rl_ez, decreasing = TRUE)
                ranked_ez[[nm]] <- rl_ez[!duplicated(names(rl_ez))]
            }
        }
    }

    if (database == "GO") {
        if (length(ranked_ens) < 2) stop("Need >= 2 contrasts with valid ENSEMBL data.")
    } else {
        if (length(ranked_ez) < 2) stop("Need >= 2 contrasts with valid ENTREZ data for KEGG.")
    }

    set.seed(42)
    cmp <- if (database == "GO") {
        clusterProfiler::compareCluster(
            geneClusters = ranked_ens,
            fun          = "gseGO",
            OrgDb        = org_db,
            keyType      = "ENSEMBL",
            ont          = go_ont,
            pvalueCutoff = 1,
            minGSSize    = min_gs,
            maxGSSize    = max_gs,
            eps          = 0,
            seed         = TRUE,
            verbose      = FALSE
        )
    } else if (database == "KEGG") {
        clusterProfiler::compareCluster(
            geneClusters = ranked_ez[lengths(ranked_ez) > 0],
            fun          = "gseKEGG",
            organism     = kegg_org,
            pvalueCutoff = 1,
            minGSSize    = min_gs,
            maxGSSize    = max_gs,
            eps          = 0,
            seed         = TRUE,
            verbose      = FALSE
        )
    } else {
        reactome_org <- if (grepl("Mm", org_db_str)) "mouse" else "human"
        clusterProfiler::compareCluster(
            geneClusters = ranked_ez[lengths(ranked_ez) > 0],
            fun          = "gsePathway",
            organism     = reactome_org,
            pvalueCutoff = 1,
            minGSSize    = min_gs,
            maxGSSize    = max_gs,
            eps          = 0,
            seed         = TRUE,
            verbose      = FALSE
        )
    }
    tryCatch(enrichplot::pairwise_termsim(cmp), error = function(e) cmp)
}
