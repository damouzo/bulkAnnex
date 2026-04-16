#!/usr/bin/env Rscript
# dge_analysis.R — DESeq2 DGE for a single contrast
# Loads pre-fitted DDS, extracts results, generates plots.

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
})

option_list <- list(
    make_option("--dds",        type = "character", help = "Path to deseq2_dds.rds"),
    make_option("--contrast",   type = "character", help = "Contrast spec: 'variable,reference,treatment'"),
    make_option("--contrast_id",type = "character", help = "Contrast identifier (used in output filenames)"),
    make_option("--vst",        type = "character", help = "Path to VST matrix TSV (for heatmaps)"),
    make_option("--lfc_thresh", type = "double",    default = 1.0,  help = "LFC threshold for significance calls and volcano [default: 1.0]"),
    make_option("--alpha",      type = "double",    default = 0.05, help = "FDR threshold [default: 0.05]"),
    make_option("--n_top",      type = "integer",   default = 50,   help = "Top N genes for heatmap [default: 50]")
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
    filter(!is.na(padj)) %>%
    mutate(
        sig        = padj < opt$alpha & abs(log2FoldChange) >= opt$lfc_thresh,
        direction  = case_when(
            padj < opt$alpha & log2FoldChange >=  opt$lfc_thresh ~ "Up",
            padj < opt$alpha & log2FoldChange <= -opt$lfc_thresh ~ "Down",
            TRUE ~ "NS"
        ),
        label      = ifelse(sig & rank(padj) <= 20, gene_name, NA_character_)
    )

dir_colours <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")

p_vol <- ggplot(vol_df, aes(x = log2FoldChange, y = -log10(padj),
                             colour = direction)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30, na.rm = TRUE) +
    geom_vline(xintercept = c(-opt$lfc_thresh, opt$lfc_thresh), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(opt$alpha), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = dir_colours) +
    labs(title  = paste0("Volcano: ", contrast_trt, " vs ", contrast_ref),
         x      = "Shrunken log2 fold change",
         y      = expression(-log[10](padj)),
         colour = NULL) +
    theme_bw(base_size = 12)

save_plot(paste0(contrast_id, "_volcano"), p_vol)

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

# ---- top genes heatmap ------------------------------------------------------
message("Plotting top genes heatmap...")

# get samples for this contrast
col_data <- as.data.frame(colData(dds))
samples_in_contrast <- rownames(col_data[col_data[[contrast_var]] %in% c(contrast_ref, contrast_trt), ])
vst_sub  <- vst_mat[, samples_in_contrast, drop = FALSE]

top_genes <- head(res_df$gene_id[!is.na(res_df$padj)], opt$n_top)
if (length(top_genes) < 2) {
    message("  Too few significant genes for heatmap — skipping.")
} else {
    hm_mat   <- vst_sub[top_genes, , drop = FALSE]
    hm_mat   <- hm_mat - rowMeans(hm_mat)  # center
    rownames(hm_mat) <- gene_name_map[rownames(hm_mat)]

    ann_col  <- data.frame(condition = col_data[samples_in_contrast, contrast_var],
                           row.names = samples_in_contrast)
    cond_lvls  <- unique(ann_col$condition)
    cond_cols  <- setNames(
        colorRampPalette(brewer.pal(min(length(cond_lvls), 9), "Set1"))(length(cond_lvls)),
        cond_lvls)
    ann_colours <- list(condition = cond_cols)

    for (ext in c("pdf", "png")) {
        f <- paste0(contrast_id, "_top_genes_heatmap.", ext)
        ht <- max(5, length(top_genes) * 0.18 + 2)
        if (ext == "png") png(f, width = 800, height = ht * 100, res = 100)
        else              pdf(f, width = 8, height = ht)
        pheatmap(hm_mat,
                 annotation_col    = ann_col,
                 annotation_colors = ann_colours,
                 cluster_rows      = TRUE,
                 cluster_cols      = TRUE,
                 show_rownames     = TRUE,
                 show_colnames     = TRUE,
                 fontsize_row      = 7,
                 fontsize_col      = 8,
                 color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                 main              = paste0("Top ", length(top_genes), " genes: ",
                                            contrast_trt, " vs ", contrast_ref))
        dev.off()
    }
}

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
pkg_versions <- c("DESeq2", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "dplyr", "ashr")
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
