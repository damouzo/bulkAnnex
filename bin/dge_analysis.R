#!/usr/bin/env Rscript
# dge_analysis.R — DESeq2 DGE for a single contrast
# Loads pre-fitted DDS, extracts results, generates plots.

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(dplyr)
    library(grid)
})

option_list <- list(
    make_option("--dds",        type = "character", help = "Path to deseq2_dds.rds"),
    make_option("--contrast",   type = "character", help = "Contrast spec: 'variable,reference,treatment'"),
    make_option("--contrast_id",type = "character", help = "Contrast identifier (used in output filenames)"),
    make_option("--vst",        type = "character", help = "Path to VST matrix TSV (for heatmaps)"),
    make_option("--lfc_thresh", type = "double",    default = 1.0,  help = "LFC threshold for significance calls and volcano [default: 1.0]"),
    make_option("--alpha",      type = "double",    default = 0.05, help = "FDR threshold [default: 0.05]"),
    make_option("--n_top",      type = "integer",   default = 50,   help = "Top N genes for the base heatmap [default: 50]"),
    make_option("--n_top_extra", type = "character", default = "100,200",
                help = "Comma-separated additional top-N heatmap sizes [default: 100,200]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- parse contrast ---------------------------------------------------------
parts <- strsplit(opt$contrast, ",")[[1]]
if (length(parts) != 3) stop("--contrast must be 'variable,reference,treatment'", call. = FALSE)
contrast_var  <- parts[1]
contrast_ref  <- parts[2]
contrast_trt  <- parts[3]
contrast_id   <- opt$contrast_id
message(sprintf("Contrast: %s  (%s: %s vs %s)", contrast_id, contrast_var, contrast_trt, contrast_ref))

# ---- load objects -----------------------------------------------------------
message("Loading DDS: ", opt$dds)
dds <- readRDS(opt$dds)

message("Loading VST matrix: ", opt$vst)
vst_raw  <- read.delim(opt$vst, check.names = FALSE, stringsAsFactors = FALSE)
gene_name_map <- setNames(vst_raw$gene_name, vst_raw$gene_id)
vst_mat  <- as.matrix(vst_raw[, !colnames(vst_raw) %in% c("gene_id", "gene_name")])
rownames(vst_mat) <- vst_raw$gene_id

# ---- extract results --------------------------------------------------------
message("Extracting DESeq2 results...")
res <- results(dds,
               contrast  = c(contrast_var, contrast_trt, contrast_ref),
               alpha     = opt$alpha,
               pAdjustMethod = "BH")

res <- lfcShrink(dds,
                 contrast = c(contrast_var, contrast_trt, contrast_ref),
                 res      = res,
                 type     = "ashr",
                 quiet    = TRUE)

res_df <- as.data.frame(res)
res_df$gene_id   <- rownames(res_df)
res_df$gene_name <- gene_name_map[res_df$gene_id]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# ---- save full results ------------------------------------------------------
out_csv <- paste0(contrast_id, "_DESeq2_results.csv")
# lfcShrink(type="ashr") may drop the 'stat' column; select only what is present
cols_wanted  <- c("gene_id", "gene_name", "baseMean", "log2FoldChange",
                  "lfcSE", "stat", "pvalue", "padj")
cols_present <- intersect(cols_wanted, colnames(res_df))
write.csv(res_df[, cols_present],
          file = out_csv, row.names = FALSE, quote = FALSE)
message("Results saved: ", out_csv)

# Preserve full in-memory precision for dashboard recomputation.
out_rds <- paste0(contrast_id, "_DESeq2_results.rds")
saveRDS(res_df[, cols_present], out_rds)
message("Results saved: ", out_rds)

# ---- summary ----------------------------------------------------------------
n_sig_up   <- sum(res_df$padj < opt$alpha & res_df$log2FoldChange > 0,  na.rm = TRUE)
n_sig_down <- sum(res_df$padj < opt$alpha & res_df$log2FoldChange < 0,  na.rm = TRUE)
message(sprintf("  Significant: %d up, %d down (padj < %.2f)", n_sig_up, n_sig_down, opt$alpha))

# ---- helper: save plot PDF + PNG -------------------------------------------
save_plot <- function(base, p, width = 8, height = 7, gg = TRUE) {
    for (ext in c("pdf", "png")) {
        f <- paste0(base, ".", ext)
        if (ext == "png") png(f, width = width * 100, height = height * 100, res = 100)
        else              pdf(f, width = width, height = height)
        if (gg) print(p) else p()
        dev.off()
    }
}

# ---- volcano plot -----------------------------------------------------------
message("Plotting volcano...")
vol_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
        sig       = padj < opt$alpha & abs(log2FoldChange) >= opt$lfc_thresh,
        padj_plot = pmax(padj, .Machine$double.xmin),
        neglog10_padj = -log10(padj_plot),
        direction = case_when(
            padj < opt$alpha & log2FoldChange >=  opt$lfc_thresh ~ "Up",
            padj < opt$alpha & log2FoldChange <= -opt$lfc_thresh ~ "Down",
            TRUE ~ "NS"
        ),
        label_base = ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, gene_id)
    )
# Label top 30 significant genes ranked by Euclidean distance from origin
# in normalised [0,1] space (equal weight to both axes regardless of scale).
label_ids_vol <- local({
    s <- vol_df[vol_df$sig, ]
    if (nrow(s) == 0) return(character(0))
    max_lfc <- max(abs(s$log2FoldChange))
    max_nlp <- max(s$neglog10_padj)
    s$score <- sqrt(
        (abs(s$log2FoldChange) / pmax(max_lfc, 1e-9))^2 +
        (s$neglog10_padj / pmax(max_nlp, 1e-9))^2
    )
    head(s$gene_id[order(s$score, decreasing = TRUE)], 30)
})
vol_df$label <- ifelse(vol_df$gene_id %in% label_ids_vol, vol_df$label_base, NA_character_)

dir_colours <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")

p_vol <- ggplot(vol_df, aes(x = log2FoldChange, y = neglog10_padj,
                             colour = direction)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(aes(label = label),
                    size = 2.5,
                    max.overlaps = Inf,
                    box.padding = 0.35,
                    point.padding = 0.2,
                    seed = 42,
                    na.rm = TRUE) +
    geom_vline(xintercept = c(-opt$lfc_thresh, opt$lfc_thresh), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(pmax(opt$alpha, .Machine$double.xmin)), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = dir_colours) +
    labs(title  = paste0("Volcano: ", contrast_trt, " vs ", contrast_ref),
         x      = "Shrunken log2 fold change",
         y      = expression(-log[10](padj)),
         colour = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.margin = margin(t = 20, r = 10, b = 5, l = 5, unit = "pt")) +
    coord_cartesian(clip = "off")

save_plot(paste0(contrast_id, "_volcano"), p_vol, width = 8, height = 8)

# ---- MA plot ----------------------------------------------------------------
message("Plotting MA...")
ma_df <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(sig   = padj < opt$alpha,
           label = ifelse(sig & rank(padj) <= 15, gene_name, NA_character_))

p_ma <- ggplot(ma_df, aes(x = log10(baseMean + 1), y = log2FoldChange,
                           colour = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    geom_hline(yintercept = 0, colour = "grey30") +
    scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey70"),
                        labels = c("TRUE" = paste0("padj < ", opt$alpha), "FALSE" = "NS")) +
    labs(title  = paste0("MA plot: ", contrast_trt, " vs ", contrast_ref),
         x      = "log10(mean normalised count + 1)",
         y      = "Shrunken log2 fold change",
         colour = NULL) +
    theme_bw(base_size = 12)

save_plot(paste0(contrast_id, "_MA"), p_ma)

# ---- top genes heatmaps (ComplexHeatmap) ------------------------------------
message("Plotting top genes heatmaps...")

# Candidate genes: significant (padj < alpha) AND |log2FC| >= lfc threshold.
# Each heatmap then takes its own top N ranked by |log2FC| (both up and down).
sig_genes_df <- res_df[
    !is.na(res_df$padj) & res_df$padj < opt$alpha &
    !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= opt$lfc_thresh,
    c("gene_id", "log2FoldChange"),
    drop = FALSE
]

select_top_genes <- function(n) {
    if (nrow(sig_genes_df) == 0) return(character(0))
    ord <- order(abs(sig_genes_df$log2FoldChange), decreasing = TRUE)
    head(sig_genes_df$gene_id[ord], n)
}

# Samples belonging to this contrast
col_data <- as.data.frame(colData(dds))
samples_in_contrast <- rownames(col_data[col_data[[contrast_var]] %in% c(contrast_ref, contrast_trt), ])
vst_sub  <- vst_mat[, samples_in_contrast, drop = FALSE]

# Column annotation (condition) - shared by all heatmaps; unsupervised throughout
ann_col <- data.frame(
    condition = factor(col_data[samples_in_contrast, contrast_var],
                       levels = unique(c(contrast_ref, contrast_trt))),
    row.names = samples_in_contrast
)
cond_lvls <- levels(ann_col$condition)
n_pal     <- max(3, min(length(cond_lvls), 9))
cond_cols <- setNames(colorRampPalette(brewer.pal(n_pal, "Set1"))(length(cond_lvls)), cond_lvls)
ha_col <- ComplexHeatmap::HeatmapAnnotation(
    condition          = ann_col$condition,
    col                = list(condition = cond_cols),
    annotation_name_gp = grid::gpar(fontsize = 8)
)

draw_top_heatmap <- function(n_top) {
    top_genes <- select_top_genes(n_top)
    if (length(top_genes) < 2) {
        message(sprintf("  Top %d: too few significant genes - skipping.", n_top))
        return(invisible(NULL))
    }
    hm_mat <- vst_sub[top_genes, , drop = FALSE]
    # Drop rows with non-finite values or zero variance (breaks clustering)
    keep   <- apply(hm_mat, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
    hm_mat <- hm_mat[keep, , drop = FALSE]
    if (nrow(hm_mat) < 2) {
        message(sprintf("  Top %d: too few variable genes - skipping.", n_top))
        return(invisible(NULL))
    }
    hm_mat <- hm_mat - rowMeans(hm_mat)  # center per row

    # Row labels: gene symbol if available, else gene_id; ensure uniqueness
    labs <- gene_name_map[rownames(hm_mat)]
    miss <- is.na(labs) | !nzchar(labs)
    labs[miss] <- rownames(hm_mat)[miss]
    rownames(hm_mat) <- make.unique(labs)

    lim <- max(stats::quantile(abs(hm_mat), 0.98, na.rm = TRUE), 1e-6)
    col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))

    # Gene names only for the small (base) heatmap; denser ones omit them.
    show_names <- n_top <= 50

    ht <- ComplexHeatmap::Heatmap(
        hm_mat,
        name                 = "Centered VST",
        col                  = col_fun,
        top_annotation       = ha_col,
        show_row_names       = show_names,
        show_column_names    = TRUE,
        row_names_gp         = grid::gpar(fontsize = 6),
        column_names_gp      = grid::gpar(fontsize = 8),
        cluster_rows         = TRUE,
        cluster_columns      = TRUE,
        column_title         = paste0("Top ", nrow(hm_mat), " genes: ",
                                      contrast_trt, " vs ", contrast_ref),
        column_title_gp      = grid::gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(title_gp  = grid::gpar(fontsize = 9),
                                    labels_gp = grid::gpar(fontsize = 8))
    )

    # Fixed 8x8 in for every heatmap (consistent png/pdf height).
    for (ext in c("pdf", "png")) {
        f <- paste0(contrast_id, "_top", n_top, "_genes_heatmap.", ext)
        if (ext == "png") png(f, width = 800, height = 800, res = 100)
        else              pdf(f, width = 8, height = 8)
        ComplexHeatmap::draw(ht)
        dev.off()
        message("  Saved: ", f)
    }
    invisible(NULL)
}

top_n_values <- unique(c(opt$n_top,
                         as.integer(trimws(strsplit(opt$n_top_extra, ",")[[1]]))))
top_n_values <- top_n_values[!is.na(top_n_values) & top_n_values > 1]
for (n_top in sort(top_n_values)) draw_top_heatmap(n_top)

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
pkg_versions <- c("DESeq2", "ggplot2", "ggrepel", "ComplexHeatmap", "circlize", "RColorBrewer", "dplyr", "ashr")
ver_lines <- c("DESEQ2_DGE:",
               paste0("    R: ", ver$R.version$version.string))
for (pkg in pkg_versions) {
    tryCatch(
        ver_lines <<- c(ver_lines, paste0("    ", pkg, ": ", as.character(packageVersion(pkg)))),
        error = function(e) NULL
    )
}
writeLines(ver_lines, "versions.yml")
message("DGE analysis complete: ", contrast_id)
