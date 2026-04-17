#!/usr/bin/env Rscript
# regenerate_qc_metrics.R — Regenerate qc_metrics.csv from raw counts + samplesheet.
# Run this if qc_metrics.csv is missing from demo_results/qc/.
# Usage: Rscript regenerate_qc_metrics.R --counts <counts.tsv> --samplesheet <ss.csv> --outdir <dir>

suppressPackageStartupMessages({
    library(DESeq2)
    library(optparse)
})

option_list <- list(
    make_option("--counts",      type = "character"),
    make_option("--samplesheet", type = "character"),
    make_option("--outdir",      type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

message("Loading counts: ", opt$counts)
counts_raw <- read.delim(opt$counts, check.names = FALSE)
gene_ids   <- counts_raw$gene_id

# Keep only numeric sample columns
sample_cols <- setdiff(colnames(counts_raw), c("gene_id", "gene_name", "tx_ids"))
counts_mat  <- round(as.matrix(counts_raw[, sample_cols, drop = FALSE]))
rownames(counts_mat) <- gene_ids

message("Loading samplesheet: ", opt$samplesheet)
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
ss <- ss[ss$sample %in% sample_cols, ]
# Reorder to match counts_mat columns
ss <- ss[match(colnames(counts_mat), ss$sample), ]

# Basic QC stats
lib_sizes  <- colSums(counts_mat)
n_detected <- colSums(counts_mat > 0)

# DESeq2 size factors
has_batch <- "batch" %in% colnames(ss) && !all(is.na(ss$batch))
design_formula <- if (has_batch) ~ batch + condition else ~ condition
coldata <- data.frame(
    condition = factor(ss$condition),
    row.names = ss$sample
)
if (has_batch) coldata$batch <- factor(ss$batch)
dds <- DESeqDataSetFromMatrix(counts_mat, colData = coldata, design = design_formula)
dds <- estimateSizeFactors(dds)
sf  <- sizeFactors(dds)

qc_metrics <- data.frame(
    sample       = ss$sample,
    condition    = ss$condition,
    library_size = lib_sizes,
    n_detected   = n_detected,
    size_factor  = round(sf, 4),
    stringsAsFactors = FALSE
)
if (has_batch) qc_metrics$batch <- ss$batch

out_path <- file.path(opt$outdir, "qc_metrics.csv")
write.csv(qc_metrics, file = out_path, row.names = FALSE, quote = FALSE)
message("Written: ", out_path)
