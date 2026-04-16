#!/usr/bin/env Rscript
# qc_analysis.R — QC plots and metrics from raw counts + DESeq2 VST

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(tidyr)
})

option_list <- list(
    make_option("--counts",      type = "character", help = "Counts matrix TSV (gene_id, gene_name, samples...)"),
    make_option("--samplesheet", type = "character", help = "Validated samplesheet CSV"),
    make_option("--prefix",      type = "character", default = "qc", help = "Output file prefix [default: qc]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- load data --------------------------------------------------------------
message("Loading counts...")
counts_raw <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts_raw) <- counts_raw$gene_id

message("Loading samplesheet...")
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
rownames(ss) <- ss$sample

# keep only samples in samplesheet, in samplesheet order
sample_cols <- ss$sample
counts_mat  <- as.matrix(counts_raw[, sample_cols])
counts_int  <- round(counts_mat)
storage.mode(counts_int) <- "integer"
gene_names  <- counts_raw[rownames(counts_int), "gene_name"]

# ---- build DESeqDataSet for QC (simple design) ------------------------------
message("Building DESeqDataSet for VST...")
has_batch <- "batch" %in% colnames(ss)
design_formula <- if (has_batch) ~ batch + condition else ~ condition

dds <- DESeqDataSetFromMatrix(
    countData = counts_int,
    colData   = ss,
    design    = design_formula
)
dds <- estimateSizeFactors(dds)
vst <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_mat <- assay(vst)

# ---- save VST matrix --------------------------------------------------------
message("Saving VST matrix...")
vst_df <- as.data.frame(vst_mat)
vst_df <- cbind(gene_id = rownames(vst_df), gene_name = gene_names, vst_df)
write.table(vst_df, file = paste0(opt$prefix, "_vst_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- QC metrics table -------------------------------------------------------
message("Computing QC metrics...")
lib_sizes   <- colSums(counts_int)
n_detected  <- colSums(counts_int > 0)
size_factors <- sizeFactors(dds)

qc_metrics <- data.frame(
    sample        = ss$sample,
    condition     = ss$condition,
    library_size  = lib_sizes,
    n_detected    = n_detected,
    size_factor   = round(size_factors, 4),
    stringsAsFactors = FALSE
)
if (has_batch) qc_metrics$batch <- ss$batch
write.csv(qc_metrics, file = paste0(opt$prefix, "_metrics.csv"), row.names = FALSE, quote = FALSE)

# ---- helper: save plot as PDF + PNG ----------------------------------------
save_plot <- function(filename_base, plot_obj = NULL, width = 8, height = 6, gg = TRUE) {
    for (ext in c("pdf", "png")) {
        f <- paste0(filename_base, ".", ext)
        if (ext == "png") png(f, width = width * 100, height = height * 100, res = 100)
        else              pdf(f, width = width, height = height)
        if (gg) print(plot_obj) else plot_obj()
        dev.off()
    }
}

# colour palette per condition
conditions     <- factor(ss$condition)
n_cond         <- nlevels(conditions)
cond_colours   <- setNames(
    colorRampPalette(brewer.pal(min(n_cond, 9), "Set1"))(n_cond),
    levels(conditions)
)

# ---- 1. Library size bar plot -----------------------------------------------
message("Plotting library sizes...")
p_lib <- ggplot(qc_metrics, aes(x = reorder(sample, library_size), y = library_size / 1e6,
                                fill = condition)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = cond_colours) +
    labs(title = "Library size per sample", x = NULL, y = "Millions of reads", fill = "Condition") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))

n_samples <- nrow(ss)
save_plot(paste0(opt$prefix, "_library_sizes"), p_lib,
          width = 8, height = max(5, n_samples * 0.3))

# ---- 2. Count distribution (log2 CPM boxplot) -------------------------------
message("Plotting count distributions...")
cpm_mat   <- sweep(counts_int, 2, lib_sizes / 1e6, FUN = "/")
log2cpm   <- log2(cpm_mat + 1)
log2cpm_long <- as.data.frame(log2cpm) %>%
    tibble::rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "log2cpm") %>%
    left_join(ss[, c("sample", "condition")], by = "sample")

p_dist <- ggplot(log2cpm_long, aes(x = sample, y = log2cpm, fill = condition)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = cond_colours) +
    labs(title = "Count distribution (log2 CPM)", x = NULL, y = "log2(CPM + 1)", fill = "Condition") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

save_plot(paste0(opt$prefix, "_count_distribution"), p_dist,
          width = max(8, n_samples * 0.4), height = 6)

# ---- 3. PCA (VST) -----------------------------------------------------------
message("Plotting PCA...")
pca_res   <- prcomp(t(vst_mat), scale. = FALSE)
pct_var   <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df    <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample    <- rownames(pca_df)
pca_df$condition <- ss[pca_df$sample, "condition"]
if (has_batch) pca_df$batch <- ss[pca_df$sample, "batch"]

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition,
                             label = sample,
                             shape = if (has_batch) batch else NULL)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
    scale_colour_manual(values = cond_colours) +
    labs(title = "PCA — VST (blind)",
         x = paste0("PC1 (", pct_var[1], "%)"),
         y = paste0("PC2 (", pct_var[2], "%)"),
         colour = "Condition", shape = if (has_batch) "Batch" else NULL) +
    theme_bw(base_size = 11)

save_plot(paste0(opt$prefix, "_pca"), p_pca, width = 8, height = 7)

# ---- 4. Sample-sample correlation heatmap -----------------------------------
message("Plotting correlation heatmap...")
cor_mat <- cor(vst_mat, method = "pearson")

ann_col <- data.frame(condition = ss$condition, row.names = ss$sample)
if (has_batch) ann_col$batch <- ss$batch
ann_colours <- list(condition = cond_colours)

for (ext in c("pdf", "png")) {
    f <- paste0(opt$prefix, "_correlation_heatmap.", ext)
    if (ext == "png") png(f, width = 900, height = 800, res = 100)
    else              pdf(f, width = 9, height = 8)
    pheatmap(cor_mat,
             annotation_col  = ann_col,
             annotation_row  = ann_col,
             annotation_colors = ann_colours,
             color           = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
             breaks          = seq(min(cor_mat), 1, length.out = 101),
             show_rownames   = TRUE,
             show_colnames   = TRUE,
             fontsize        = 8,
             main            = "Sample-sample Pearson correlation (VST)")
    dev.off()
}

# ---- 5. Dispersion estimate plot -------------------------------------------
message("Fitting DESeq2 dispersion for QC plot...")
dds_disp <- estimateDispersions(dds, quiet = TRUE)

for (ext in c("pdf", "png")) {
    f <- paste0(opt$prefix, "_dispersion.", ext)
    if (ext == "png") png(f, width = 700, height = 600, res = 100)
    else              pdf(f, width = 7, height = 6)
    plotDispEsts(dds_disp, main = "DESeq2 dispersion estimates")
    dev.off()
}

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
pkg_versions <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "dplyr", "tidyr")
ver_lines <- c("DESEQ2_QC:",
               paste0("    R: ", ver$R.version$version.string))
for (pkg in pkg_versions) {
    ver_lines <- c(ver_lines, paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
}
writeLines(ver_lines, "versions.yml")
message("QC analysis complete.")
