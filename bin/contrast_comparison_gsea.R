#!/usr/bin/env Rscript
# contrast_comparison_gsea.R — Cross-contrast GSEA pathway comparison
# Generates, for one or more databases, Venn diagrams, UpSet plots and a NES
# heatmap over the pathways that are significantly enriched across contrasts,
# plus a summary CSV of the pathways common to ALL contrasts.
#
# Input: *_GO_all_gsea.csv, *_KEGG_gsea.csv, *_Reactome_gsea.csv staged in the
#        Nextflow work directory.
# Output: contrast_comparison_gsea_<TAG>_{venn,upset,heatmap}.{pdf,png},
#         contrast_comparison_gsea_common_pathways_<TAG>.csv, versions.yml

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(ggVennDiagram)
    library(ComplexUpset)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(dplyr)
})

option_list <- list(
    make_option("--padj", type = "double", default = 0.05,
                help = "Adjusted p-value threshold for enriched pathways [default: %default]"),
    make_option("--databases", type = "character", default = "GO:BP",
                help = "Comma-separated: GO:BP, GO:MF, GO:CC, KEGG, Reactome [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Discover staged GSEA CSV files -----------------------------------------
go_files       <- sort(Sys.glob("*_GO_all_gsea.csv"))
kegg_files     <- sort(Sys.glob("*_KEGG_gsea.csv"))
reactome_files <- sort(Sys.glob("*_Reactome_gsea.csv"))

contrasts_of <- function(files, suffix) {
    if (length(files) == 0) return(character(0))
    sub(paste0(suffix, "$"), "", basename(files))
}
go_contrasts       <- contrasts_of(go_files,       "_GO_all_gsea.csv")
kegg_contrasts     <- contrasts_of(kegg_files,     "_KEGG_gsea.csv")
reactome_contrasts <- contrasts_of(reactome_files, "_Reactome_gsea.csv")

all_contrasts <- unique(c(go_contrasts, kegg_contrasts, reactome_contrasts))

write_skip_versions <- function(msg) {
    writeLines(c("CONTRAST_COMPARISON_GSEA:",
                 paste0("    status: skipped (", msg, ")")), "versions.yml")
}

if (length(all_contrasts) < 2) {
    message("Contrast comparison (GSEA) requires >= 2 contrasts — skipping.")
    write_skip_versions(paste0(length(all_contrasts), " contrast(s)"))
    quit(save = "no", status = 0)
}

message(sprintf("Cross-contrast GSEA comparison: %d contrasts (padj < %g)",
                length(all_contrasts), opt$padj))
for (cid in all_contrasts) message("  ", cid)

# ---- Helpers ----------------------------------------------------------------
read_one <- function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"padj" %in% colnames(df) && "p.adjust" %in% colnames(df)) df$padj <- df$p.adjust
    if (!"padj" %in% colnames(df) && "pval"     %in% colnames(df)) df$padj <- p.adjust(df$pval, method = "BH")
    if (!"ID"   %in% colnames(df) && "Description" %in% colnames(df)) df$ID <- df$Description
    if (!"NES"  %in% colnames(df)) df$NES <- NA_real_
    df
}

load_table <- function(type, contrast) {
    files  <- switch(type, GO = go_files,       KEGG = kegg_files,     Reactome = reactome_files)
    contr  <- switch(type, GO = go_contrasts,   KEGG = kegg_contrasts, Reactome = reactome_contrasts)
    idx    <- match(contrast, contr)
    if (is.na(idx)) return(NULL)
    read_one(files[idx])
}

parse_db <- function(spec) {
    spec <- trimws(spec)
    if (toupper(spec) == "KEGG")    return(list(type = "KEGG",     ont = NA_character_, tag = "KEGG",     label = "KEGG"))
    if (toupper(spec) == "REACTOME") return(list(type = "Reactome", ont = NA_character_, tag = "Reactome", label = "Reactome"))
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
    ont   <- if (length(parts) >= 2 && nzchar(parts[2])) toupper(parts[2]) else "BP"
    list(type = "GO", ont = ont, tag = paste0("GO_", ont), label = paste0("GO ", ont))
}

db_specs <- lapply(strsplit(opt$databases, ",")[[1]], parse_db)

save_plot <- function(p, prefix, width = 10, height = 8) {
    pdf_file <- paste0(prefix, ".pdf")
    png_file <- paste0(prefix, ".png")
    pdf(pdf_file, width = width, height = height)
    if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) draw(p) else print(p)
    dev.off()
    png(png_file, width = width, height = height, units = "in", res = 150)
    if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) draw(p) else print(p)
    dev.off()
    message("Saved: ", pdf_file, " and ", png_file)
}

make_placeholder <- function(msg) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, colour = "grey40") +
        theme_void()
}

trunc_lab <- function(x, n = 45) {
    ifelse(is.na(x), "", ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x))
}

# ---- Per-database comparison ------------------------------------------------
for (spec in db_specs) {
    message("\n=== Database: ", spec$label, " ===")

    dfs <- setNames(lapply(all_contrasts, function(cid) {
        df <- load_table(spec$type, cid)
        if (is.null(df)) return(NULL)
        if (spec$type == "GO" && "ontology" %in% colnames(df)) {
            df <- df[toupper(df$ontology) == spec$ont, , drop = FALSE]
        }
        df[!is.na(df$padj) & df$padj < opt$padj, , drop = FALSE]
    }), all_contrasts)
    dfs <- dfs[!sapply(dfs, is.null)]

    if (length(dfs) < 2) {
        message("Skipping ", spec$tag, ": < 2 contrasts with data.")
        next
    }

    sets <- lapply(dfs, function(d) if (nrow(d) == 0) character(0) else unique(d$ID))

    labels <- local({
        all_df <- bind_rows(lapply(dfs, function(d) d[, c("ID", "Description"), drop = FALSE]))
        if (nrow(all_df) == 0) return(character(0))
        all_df <- all_df[!duplicated(all_df$ID), , drop = FALSE]
        setNames(all_df$Description, all_df$ID)
    })

    n_sig <- sapply(sets, length)
    message("Significant pathways per contrast: ",
            paste(paste0(names(n_sig), "=", n_sig), collapse = ", "))

    # --- Venn ---------------------------------------------------------------
    venn_p <- if (length(sets) > 5) {
        make_placeholder(paste0("Venn Diagram skipped:\n", length(sets),
                                " contrasts exceed the 5-set limit.\nSee UpSet plot instead."))
    } else if (all(lengths(sets) == 0)) {
        make_placeholder("No significant pathways found.")
    } else {
        ggVennDiagram(sets, label = "both", set_color = rep("#4d9221", length(sets))) +
            scale_fill_gradient(low = "white", high = "#4d9221") +
            scale_color_manual(values = rep("#4d9221", length(sets))) +
            ggtitle(paste0("Common GSEA pathways — ", spec$label)) +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    }
    save_plot(venn_p, paste0("contrast_comparison_gsea_", spec$tag, "_venn"),
              width = 8, height = 7)

    # --- UpSet --------------------------------------------------------------
    upset_p <- local({
        all_paths <- unique(unlist(sets))
        if (length(all_paths) == 0) return(make_placeholder("No significant pathways found."))
        mat <- as.data.frame(
            sapply(names(sets), function(cid) all_paths %in% sets[[cid]]),
            stringsAsFactors = FALSE)
        mat[] <- lapply(mat, as.logical)
        upset(mat, intersect = names(sets), name = "Contrast",
              base_annotations = list(
                  "Intersection\nsize" = intersection_size(
                      counts  = TRUE,
                      mapping = aes(fill = !!sym(names(sets)[1]))) +
                      scale_fill_manual(values = c("TRUE" = "#4d9221", "FALSE" = "#4d9221"),
                                        guide = "none")),
              set_sizes = upset_set_size(geom = geom_bar(fill = "#4d9221")),
              width_ratio = 0.2, sort_sets = FALSE,
              sort_intersections_by = "cardinality", n_intersections = 40) +
            ggtitle(paste0("Common GSEA pathways — ", spec$label)) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
    })
    save_plot(upset_p, paste0("contrast_comparison_gsea_", spec$tag, "_upset"),
              width = 12, height = 7)

    # --- NES heatmap (pathways common to ALL contrasts) ---------------------
    common <- Reduce(intersect, sets)
    message("Common pathways (all contrasts): ", length(common))

    heatmap_p <- local({
        if (length(common) < 2) {
            return(make_placeholder(paste0(
                "Not enough pathways enriched in all contrasts\n(n = ", length(common), ").")))
        }
        nes_mat <- sapply(names(dfs), function(cid) {
            d <- dfs[[cid]]
            vec <- setNames(d$NES, d$ID)
            vec[common]
        })
        nes_mat <- as.matrix(nes_mat)
        rownames(nes_mat) <- common

        labs <- labels[common]
        miss <- is.na(labs) | !nzchar(labs)
        labs[miss] <- common[miss]
        rownames(nes_mat) <- make.unique(trunc_lab(labs, 45))

        lim     <- max(stats::quantile(abs(nes_mat), 0.98, na.rm = TRUE), 1e-6)
        col_fun <- colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))

        Heatmap(
            nes_mat,
            name              = "NES",
            col               = col_fun,
            na_col            = "#cccccc",
            show_row_names    = nrow(nes_mat) <= 50,
            row_names_gp      = gpar(fontsize = 7),
            row_names_max_width = unit(7, "cm"),
            column_names_gp   = gpar(fontsize = 10),
            cluster_rows      = TRUE,
            cluster_columns   = FALSE,
            column_title      = paste0("Common GSEA pathways — ", spec$label,
                                       "\n(n = ", nrow(nes_mat), " x ", ncol(nes_mat), " contrasts)"),
            column_title_gp   = gpar(fontsize = 13, fontface = "bold"),
            heatmap_legend_param = list(title_gp  = gpar(fontsize = 10),
                                        labels_gp = gpar(fontsize = 9))
        )
    })
    save_plot(heatmap_p, paste0("contrast_comparison_gsea_", spec$tag, "_heatmap"),
              width = 9, height = 8)

    # --- Common pathways summary CSV ----------------------------------------
    csv_file <- paste0("contrast_comparison_gsea_common_pathways_", spec$tag, ".csv")
    if (length(common) > 0) {
        nes_cols <- sapply(names(dfs), function(cid) {
            round(setNames(dfs[[cid]]$NES, dfs[[cid]]$ID)[common], 4)
        })
        padj_cols <- sapply(names(dfs), function(cid) {
            signif(setNames(dfs[[cid]]$padj, dfs[[cid]]$ID)[common], 3)
        })
        if (is.null(dim(nes_cols))) {
            nes_cols  <- matrix(nes_cols,  nrow = 1, dimnames = list(common, names(dfs)))
            padj_cols <- matrix(padj_cols, nrow = 1, dimnames = list(common, names(dfs)))
        }
        colnames(nes_cols)  <- paste0("NES_",  names(dfs))
        colnames(padj_cols) <- paste0("padj_", names(dfs))
        summary_df <- data.frame(
            pathway_id  = common,
            description = labels[common],
            nes_cols,
            padj_cols,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        write.csv(summary_df, csv_file, row.names = FALSE, quote = FALSE)
    } else {
        write.csv(data.frame(pathway_id = character(), description = character()),
                  csv_file, row.names = FALSE, quote = FALSE)
    }
    message("Saved: ", csv_file)
}

# ---- versions.yml -----------------------------------------------------------
versions_content <- c(
    "CONTRAST_COMPARISON_GSEA:",
    paste0("    R: ", R.version$major, ".", R.version$minor),
    paste0("    ggVennDiagram: ", packageVersion("ggVennDiagram")),
    paste0("    ComplexUpset:  ", packageVersion("ComplexUpset")),
    paste0("    ComplexHeatmap: ", packageVersion("ComplexHeatmap")),
    paste0("    databases: ", opt$databases),
    paste0("    n_contrasts: ", length(all_contrasts)),
    paste0("    padj_threshold: ", opt$padj)
)
writeLines(versions_content, "versions.yml")
message("\ncontrast_comparison_gsea.R completed successfully.")
