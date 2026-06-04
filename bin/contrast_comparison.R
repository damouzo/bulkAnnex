#!/usr/bin/env Rscript
# contrast_comparison.R — Cross-contrast DEG comparison
# Generates Venn diagrams, UpSet plots, and Log2FC Spearman correlation
# heatmaps for UP- and DOWN-regulated DEGs across all contrasts.
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
})

option_list <- list(
    make_option("--padj", type = "double",  default = 0.05,
                help = "Adjusted p-value threshold for DEG calls [default: %default]"),
    make_option("--lfc",  type = "double",  default = 0.5,
                help = "Absolute log2FC threshold for DEG calls [default: %default]")
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
                mapping   = aes(fill = !!sym(names(sets)[1]))   # dummy fill overridden below
            ) +
            scale_fill_manual(values = colour, guide = "none")
        ),
        set_sizes = upset_set_size(
            geom = geom_bar(fill = colour)
        ),
        width_ratio   = 0.2,
        sort_sets     = FALSE,
        sort_intersections_by = "degree"
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

# ============================================================
# SECTION 4: Common DEGs summary CSV
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
        lfc_cols <- sapply(names(dge_list), function(cid) {
            df  <- dge_list[[cid]]
            vec <- setNames(df$log2FoldChange, df$gene_id)
            round(vec[gene_ids], 4)
        })
        padj_cols <- sapply(names(dge_list), function(cid) {
            df  <- dge_list[[cid]]
            vec <- setNames(df$padj, df$gene_id)
            signif(vec[gene_ids], 3)
        })
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
# SECTION 5: versions.yml
# ============================================================
versions_content <- c(
    "CONTRAST_COMPARISON:",
    paste0("    R: ", R.version$major, ".", R.version$minor),
    paste0("    ggVennDiagram: ", packageVersion("ggVennDiagram")),
    paste0("    ComplexUpset: ",  packageVersion("ComplexUpset")),
    paste0("    ComplexHeatmap: ", packageVersion("ComplexHeatmap")),
    paste0("    n_contrasts: ", length(dge_files)),
    paste0("    padj_threshold: ", opt$padj),
    paste0("    lfc_threshold: ", opt$lfc)
)
writeLines(versions_content, "versions.yml")
message("\ncontrast_comparison.R completed successfully.")
