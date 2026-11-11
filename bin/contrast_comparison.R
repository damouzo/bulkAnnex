#!/usr/bin/env Rscript
# contrast_comparison.R — Cross-contrast DEG comparison
# Generates Venn diagrams, UpSet plots, and Log2FC Spearman correlation
# heatmaps for UP- and DOWN-regulated DEGs across all contrasts, plus a
# gene heatmap (log2FC) and a GO BP ORA dotplot over the genes that are DEG
# in ALL contrasts (either direction).
#
# Input: all *_DESeq2_results.csv files staged in the Nextflow work directory.
# Output: contrast_comparison_*.{pdf,png}, contrast_comparison_common_degs.csv,
#         versions.yml

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(ggVennDiagram)
    library(ComplexUpset)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(dplyr)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
})

option_list <- list(
    make_option("--padj", type = "double",  default = 0.05,
                help = "Adjusted p-value threshold for DEG calls [default: %default]"),
    make_option("--lfc",  type = "double",  default = 1.0,
                help = "Absolute log2FC threshold for DEG calls [default: %default]"),
    make_option("--organism", type = "character", default = "auto",
                help = "Organism: human, mouse, or auto (detect from gene IDs) [default: %default]"),
    make_option("--ora_topn", type = "integer", default = 10,
                help = "Number of top GO BP terms in the ORA dotplot [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Discover DGE CSV files staged in the work directory --------------------
dge_files <- sort(Sys.glob("*_DESeq2_results.csv"))

if (length(dge_files) < 2) {
    msg <- if (length(dge_files) == 0) "No *_DESeq2_results.csv files found."
           else "Only 1 contrast — contrast comparison requires >= 2 contrasts."
    message(msg, " Skipping contrast_comparison.")
    writeLines(c("CONTRAST_COMPARISON:",
                 paste0("    status: skipped (", length(dge_files), " contrast(s))")),
               "versions.yml")
    quit(save = "no", status = 0)
}

contrast_ids <- sub("_DESeq2_results\\.csv$", "", basename(dge_files))
message(sprintf("Cross-contrast comparison: %d contrasts (padj < %g, |lfc| > %g)",
                length(dge_files), opt$padj, opt$lfc))
for (cid in contrast_ids) message("  ", cid)

# ---- Read DGE results -------------------------------------------------------
dge_list <- setNames(
    lapply(dge_files, function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
        # Normalise column names — some runs may use 'log2FC' instead of 'log2FoldChange'
        if (!"log2FoldChange" %in% colnames(df) && "log2FC" %in% colnames(df)) {
            df$log2FoldChange <- df$log2FC
        }
        df
    }),
    contrast_ids
)

# ---- Build directional gene sets --------------------------------------------
get_gene_sets <- function(direction) {
    lapply(setNames(contrast_ids, contrast_ids), function(cid) {
        df        <- dge_list[[cid]]
        keep_padj <- !is.na(df$padj) & df$padj < opt$padj
        keep_lfc  <- if (direction == "up") df$log2FoldChange >  opt$lfc
                     else                   df$log2FoldChange < -opt$lfc
        df$gene_id[keep_padj & keep_lfc]
    })
}

sets_up   <- get_gene_sets("up")
sets_down <- get_gene_sets("down")

n_up   <- sapply(sets_up,   length)
n_down <- sapply(sets_down, length)
message("UP DEGs per contrast:   ", paste(paste0(names(n_up),   "=", n_up),   collapse = ", "))
message("DOWN DEGs per contrast: ", paste(paste0(names(n_down), "=", n_down), collapse = ", "))

# ---- Helper: save plot in PDF + PNG -----------------------------------------
save_plot <- function(p, prefix, width = 10, height = 8) {
    pdf_file <- paste0(prefix, ".pdf")
    png_file <- paste0(prefix, ".png")

    # PDF
    pdf(pdf_file, width = width, height = height)
    if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) {
        draw(p)
    } else {
        print(p)
    }
    dev.off()

    # PNG
    png(png_file, width = width, height = height, units = "in", res = 150)
    if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) {
        draw(p)
    } else {
        print(p)
    }
    dev.off()

    message("Saved: ", pdf_file, " and ", png_file)
}

# ---- Helper: placeholder when plot cannot be drawn --------------------------
make_placeholder <- function(msg) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, colour = "grey40") +
        theme_void()
}

# ============================================================
# SECTION 1: Venn Diagrams
# ============================================================
make_venn <- function(sets, direction) {
    n_sets <- length(sets)
    label  <- if (direction == "up") "UP-regulated DEGs" else "DOWN-regulated DEGs"
    colour <- if (direction == "up") "#d73027" else "#4575b4"

    # Venn diagrams are only readable up to 5 sets.
    if (n_sets > 5) {
        message("Skipping Venn for ", direction, ": ", n_sets, " contrasts > 5 (use UpSet instead).")
        return(make_placeholder(
            paste0("Venn Diagram skipped:\n", n_sets,
                   " contrasts exceed the 5-set limit.\nSee UpSet plot instead.")))
    }

    all_empty <- all(sapply(sets, length) == 0)
    if (all_empty) {
        return(make_placeholder(paste0("No significant ", label, " found.")))
    }

    ggVennDiagram(sets, label = "both",
                  set_color = rep(colour, n_sets)) +
        scale_fill_gradient(low = "white", high = colour) +
        scale_color_manual(values = rep(colour, n_sets)) +
        ggtitle(label) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
}

message("\n--- Venn Diagrams ---")
venn_up   <- make_venn(sets_up,   "up")
venn_down <- make_venn(sets_down, "down")
save_plot(venn_up,   "contrast_comparison_venn_up",   width = 8, height = 7)
save_plot(venn_down, "contrast_comparison_venn_down", width = 8, height = 7)

# ============================================================
# SECTION 2: UpSet Plots
# ============================================================
make_upset <- function(sets, direction) {
    label  <- if (direction == "up") "UP-regulated DEGs" else "DOWN-regulated DEGs"
    colour <- if (direction == "up") "#d73027" else "#4575b4"

    all_genes <- unique(unlist(sets))
    if (length(all_genes) == 0) {
        return(make_placeholder(paste0("No significant ", label, " found.")))
    }

    # Build binary membership data frame (rows = genes, cols = contrasts)
    mat <- as.data.frame(
        sapply(names(sets), function(cid) all_genes %in% sets[[cid]]),
        stringsAsFactors = FALSE
    )
    rownames(mat) <- all_genes

    # ComplexUpset requires columns to be logical
    mat[] <- lapply(mat, as.logical)

    upset(
        mat,
        intersect = names(sets),
        name      = "Contrast",
        base_annotations = list(
            "Intersection\nsize" = intersection_size(
                counts    = TRUE,
                mapping   = aes(fill = !!sym(names(sets)[1]))
            ) +
            # Provide values for both logical levels (TRUE/FALSE) — a single value
            # causes "Insufficient values in manual scale" with ggplot2 >= 3.4.
            scale_fill_manual(
                values = c("TRUE" = colour, "FALSE" = colour),
                guide  = "none"
            )
        ),
        set_sizes = upset_set_size(
            geom = geom_bar(fill = colour)
        ),
        width_ratio   = 0.2,
        sort_sets     = FALSE,
        # "cardinality" shows the largest intersections first and avoids
        # rendering all 2^n combinations when n is large (e.g. 10 contrasts).
        sort_intersections_by = "cardinality",
        n_intersections       = 40   # cap at 40 bars to keep the plot readable
    ) +
    ggtitle(label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
}

message("\n--- UpSet Plots ---")
upset_up   <- make_upset(sets_up,   "up")
upset_down <- make_upset(sets_down, "down")
save_plot(upset_up,   "contrast_comparison_upset_up",   width = 12, height = 7)
save_plot(upset_down, "contrast_comparison_upset_down", width = 12, height = 7)

# ============================================================
# SECTION 3: Log2FC Spearman Correlation Heatmaps
# ============================================================
make_lfc_heatmap <- function(sets, dge_full, direction) {
    label <- if (direction == "up") "UP-regulated" else "DOWN-regulated"

    # Union of significant genes for this direction
    sig_genes <- unique(unlist(sets))

    if (length(sig_genes) < 2) {
        message("Skipping LFC heatmap for ", direction,
                ": insufficient significant genes (n=", length(sig_genes), ")")
        return(NULL)
    }

    # Build LFC matrix: rows = sig_genes, cols = contrasts
    # Use the full (not filtered) LFC for each gene in each contrast.
    lfc_mat <- sapply(names(dge_full), function(cid) {
        df  <- dge_full[[cid]]
        vec <- setNames(df$log2FoldChange, df$gene_id)
        vec[sig_genes]   # NA where gene absent in this contrast
    })
    rownames(lfc_mat) <- sig_genes

    # Retain only genes with non-NA LFC in at least 2 contrasts
    n_present <- rowSums(!is.na(lfc_mat))
    lfc_mat   <- lfc_mat[n_present >= 2, , drop = FALSE]

    if (nrow(lfc_mat) < 2 || ncol(lfc_mat) < 2) {
        message("Skipping LFC heatmap for ", direction, ": insufficient overlap.")
        return(NULL)
    }

    # Spearman correlation between contrasts
    cor_mat <- cor(lfc_mat, method = "spearman", use = "pairwise.complete.obs")

    # Replace NaN (single non-NA obs) with NA
    cor_mat[is.nan(cor_mat)] <- NA

    col_fun <- colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027"))

    ht <- Heatmap(
        cor_mat,
        name              = "Spearman r",
        col               = col_fun,
        column_title      = paste0(label, " DEGs — Spearman Log2FC correlation\n",
                                   "(n = ", nrow(lfc_mat), " genes)"),
        column_title_gp   = gpar(fontsize = 13, fontface = "bold"),
        row_names_gp      = gpar(fontsize = 10),
        column_names_gp   = gpar(fontsize = 10),
        show_row_names    = TRUE,
        show_column_names = TRUE,
        cluster_rows      = TRUE,
        cluster_columns   = TRUE,
        rect_gp           = gpar(col = "white", lwd = 1),
        na_col            = "#cccccc",
        cell_fun          = function(j, i, x, y, width, height, fill) {
            if (!is.na(cor_mat[i, j])) {
                grid.text(sprintf("%.2f", cor_mat[i, j]), x, y,
                          gp = gpar(fontsize = 8, col = "black"))
            }
        }
    )
    ht
}

message("\n--- Log2FC Correlation Heatmaps ---")
ht_up   <- make_lfc_heatmap(sets_up,   dge_list, "up")
ht_down <- make_lfc_heatmap(sets_down, dge_list, "down")

if (!is.null(ht_up)) {
    save_plot(ht_up,   "contrast_comparison_lfc_heatmap_up",   width = 9, height = 8)
}
if (!is.null(ht_down)) {
    save_plot(ht_down, "contrast_comparison_lfc_heatmap_down", width = 9, height = 8)
}

# ---- Common DEGs across ALL contrasts (significant in either direction) -----
# Single set used by both the gene heatmap and the ORA dotplot: genes that are
# DEG (padj < cutoff AND |log2FC| >= cutoff, up or down) in EVERY contrast.
any_sets <- lapply(setNames(contrast_ids, contrast_ids), function(cid) {
    unique(c(sets_up[[cid]], sets_down[[cid]]))
})
common_any <- if (length(any_sets) >= 2) Reduce(intersect, any_sets) else character(0)
message("Common DEGs across all contrasts (any direction): ", length(common_any))

# gene_id -> gene_name map (first occurrence across contrasts)
gene_name_map_all <- local({
    all_df <- bind_rows(lapply(dge_list, function(df) df[, c("gene_id", "gene_name"), drop = FALSE]))
    all_df <- all_df[!duplicated(all_df$gene_id), ]
    setNames(all_df$gene_name, all_df$gene_id)
})

# ============================================================
# SECTION 4: Common DEGs Gene Heatmap (log2FC across all contrasts)
# ============================================================
make_common_heatmap <- function(genes) {
    if (length(genes) < 2) {
        return(make_placeholder(paste0(
            "Not enough common DEGs across all contrasts\n(n = ", length(genes), ").")))
    }

    # log2FC matrix: rows = common genes, cols = all contrasts
    lfc_mat <- sapply(names(dge_list), function(cid) {
        df  <- dge_list[[cid]]
        vec <- setNames(df$log2FoldChange, df$gene_id)
        vec[genes]
    })
    rownames(lfc_mat) <- genes

    # Row labels: gene symbol if available, else gene_id
    labs <- gene_name_map_all[genes]
    miss <- is.na(labs) | !nzchar(labs)
    labs[miss] <- genes[miss]
    rownames(lfc_mat) <- make.unique(labs)

    lim     <- max(stats::quantile(abs(lfc_mat), 0.98, na.rm = TRUE), 1e-6)
    col_fun <- colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))

    Heatmap(
        lfc_mat,
        name              = "log2FC",
        col               = col_fun,
        na_col            = "#cccccc",
        show_row_names    = nrow(lfc_mat) <= 50,
        row_names_gp      = gpar(fontsize = 7),
        column_names_gp   = gpar(fontsize = 9),
        cluster_rows      = TRUE,
        cluster_columns   = FALSE,
        column_title      = paste0("Common DEGs across all contrasts (n = ",
                                   nrow(lfc_mat), " genes x ", ncol(lfc_mat), " contrasts)"),
        column_title_gp   = gpar(fontsize = 13, fontface = "bold"),
        heatmap_legend_param = list(title_gp  = gpar(fontsize = 10),
                                    labels_gp = gpar(fontsize = 9))
    )
}

message("\n--- Common DEGs Gene Heatmap ---")
ht_common <- make_common_heatmap(common_any)
save_plot(ht_common, "contrast_comparison_gene_heatmap_common", width = 8, height = 8)

# ============================================================
# SECTION 5: Common DEGs ORA (GO BP)
# ============================================================
detect_organism <- function(gene_ids) {
    ids <- gene_ids[!is.na(gene_ids)]
    if (length(ids) > 0 && any(grepl("^ENSMUSG", ids))) "org.Mm.eg.db" else "org.Hs.eg.db"
}

ora_org_str <- if (!is.null(opt$organism) && opt$organism %in% c("human", "mouse")) {
    if (opt$organism == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
} else {
    detect_organism(common_any)
}
message("ORA organism: ", ora_org_str)

make_common_ora <- function(genes, top_n = 10) {
    if (length(genes) < 2) {
        return(make_placeholder(paste0(
            "Not enough common DEGs across all contrasts\n(n = ", length(genes), ").")))
    }
    ego <- tryCatch(
        enrichGO(gene          = genes,
                 OrgDb         = get(ora_org_str),
                 keyType       = "ENSEMBL",
                 ont           = "BP",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1,
                 pAdjustMethod = "BH",
                 minGSSize     = 5,
                 maxGSSize     = 500,
                 readable      = TRUE),
        error = function(e) { message("enrichGO: ", e$message); NULL })

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        return(make_placeholder(paste0(
            "No GO BP terms enriched for the ", length(genes), " common DEGs.")))
    }

    enrichplot::dotplot(ego,
                        showCategory = min(top_n, nrow(as.data.frame(ego))),
                        title        = paste0("ORA GO BP - ", length(genes),
                                              " common DEGs (all contrasts)")) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

message("\n--- Common DEGs ORA (GO BP) ---")
p_ora_common <- make_common_ora(common_any, top_n = opt$ora_topn)
save_plot(p_ora_common, "contrast_comparison_ora_go_bp_common", width = 9, height = 7)

# ============================================================
# SECTION 6: Common DEGs summary CSV
# ============================================================
message("\n--- Common DEGs summary ---")

# Genes significant in ALL contrasts (strict intersection)
common_up   <- Reduce(intersect, sets_up[sapply(sets_up,   length) > 0])
common_down <- Reduce(intersect, sets_down[sapply(sets_down, length) > 0])

if ((length(common_up) + length(common_down)) > 0) {
    # Annotate with gene_name from the first contrast that has each gene
    get_gene_name <- function(gene_ids) {
        all_df <- bind_rows(lapply(dge_list, function(df) {
            df[, c("gene_id", "gene_name"), drop = FALSE]
        }))
        all_df <- all_df[!duplicated(all_df$gene_id), ]
        name_map <- setNames(all_df$gene_name, all_df$gene_id)
        name_map[gene_ids]
    }

    make_summary_df <- function(gene_ids, direction) {
        if (length(gene_ids) == 0) return(data.frame())
        # Use cbind, NOT sapply: with a single gene, sapply collapses the
        # per-contrast values into a 1-D vector, so colnames<- below fails with
        # "attempt to set 'colnames' on an object with less than two dimensions".
        # cbind always yields a genes x contrasts matrix.
        lfc_cols <- do.call(cbind, lapply(names(dge_list), function(cid) {
            df  <- dge_list[[cid]]
            vec <- setNames(df$log2FoldChange, df$gene_id)
            round(vec[gene_ids], 4)
        }))
        padj_cols <- do.call(cbind, lapply(names(dge_list), function(cid) {
            df  <- dge_list[[cid]]
            vec <- setNames(df$padj, df$gene_id)
            signif(vec[gene_ids], 3)
        }))
        colnames(lfc_cols)  <- paste0("lfc_",  names(dge_list))
        colnames(padj_cols) <- paste0("padj_", names(dge_list))
        data.frame(
            gene_id   = gene_ids,
            gene_name = get_gene_name(gene_ids),
            direction = direction,
            lfc_cols,
            padj_cols,
            stringsAsFactors = FALSE
        )
    }

    summary_df <- bind_rows(
        make_summary_df(common_up,   "up"),
        make_summary_df(common_down, "down")
    )

    write.csv(summary_df, "contrast_comparison_common_degs.csv",
              row.names = FALSE, quote = FALSE)
    message("Common DEGs (all contrasts): ",
            length(common_up), " UP, ", length(common_down), " DOWN")
    message("Saved: contrast_comparison_common_degs.csv")
} else {
    message("No genes in common across ALL contrasts.")
    # Write empty file so publishDir pattern is satisfied
    write.csv(data.frame(gene_id = character(), gene_name = character(),
                         direction = character()),
              "contrast_comparison_common_degs.csv",
              row.names = FALSE, quote = FALSE)
}

# ============================================================
# SECTION 7: versions.yml
# ============================================================
versions_content <- c(
    "CONTRAST_COMPARISON:",
    paste0("    R: ", R.version$major, ".", R.version$minor),
    paste0("    ggVennDiagram: ", packageVersion("ggVennDiagram")),
    paste0("    ComplexUpset:  ", packageVersion("ComplexUpset")),
    paste0("    ComplexHeatmap: ", packageVersion("ComplexHeatmap")),
    paste0("    clusterProfiler: ", packageVersion("clusterProfiler")),
    paste0("    enrichplot: ", packageVersion("enrichplot")),
    paste0("    org_db: ", ora_org_str),
    paste0("    n_contrasts: ", length(dge_files)),
    paste0("    n_common_genes: ", length(common_any)),
    paste0("    ora_topn: ", opt$ora_topn),
    paste0("    padj_threshold: ", opt$padj),
    paste0("    lfc_threshold: ", opt$lfc)
)
writeLines(versions_content, "versions.yml")
message("\ncontrast_comparison.R completed successfully.")
