#!/usr/bin/env Rscript
# normalization.R — DESeq2 VST normalization and full model fitting
# Saves the fitted DESeqDataSet (.rds) for downstream DGE reuse.

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
})

option_list <- list(
    make_option("--counts",      type = "character", help = "Counts matrix TSV (gene_id, gene_name, samples...)"),
    make_option("--samplesheet", type = "character", help = "Validated samplesheet CSV"),
    make_option("--prefix",      type = "character", default = "deseq2",
                help = "Output file prefix [default: deseq2]"),
    make_option("--min_counts",  type = "integer",   default = 10,
                help = "Min count threshold for gene retention [default: 10]"),
    make_option("--min_samples", type = "integer",   default = 2,
                help = "Min samples with >= min_counts to keep a gene [default: 2]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- load data --------------------------------------------------------------
message("Loading counts...")
counts_raw <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts_raw) <- counts_raw$gene_id

message("Loading samplesheet...")
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
rownames(ss) <- ss$sample

# enforce factor levels: reference level for condition = first alphabetical
ss$condition <- factor(ss$condition)

has_batch <- "batch" %in% colnames(ss)
if (has_batch) {
    ss$batch <- factor(ss$batch)
    message("Batch column detected. Design: ~ batch + condition")
    design_formula <- ~ batch + condition
} else {
    message("No batch column. Design: ~ condition")
    design_formula <- ~ condition
}

# ---- build counts matrix ----------------------------------------------------
sample_cols <- ss$sample
counts_mat  <- as.matrix(counts_raw[, sample_cols])
counts_int  <- round(counts_mat)
storage.mode(counts_int) <- "integer"

# filter: keep genes with >= min_counts in at least min_samples samples
keep <- rowSums(counts_int >= opt$min_counts) >= opt$min_samples
message(sprintf("Filtering: %d / %d genes pass (>= %d counts in >= %d samples).",
                sum(keep), nrow(counts_int), opt$min_counts, opt$min_samples))
counts_int <- counts_int[keep, ]

# store gene_name as row metadata
gene_meta <- data.frame(
    gene_id   = rownames(counts_int),
    gene_name = counts_raw[rownames(counts_int), "gene_name"],
    stringsAsFactors = FALSE
)

# ---- DESeq2 -----------------------------------------------------------------
message("Building DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
    countData = counts_int,
    colData   = ss,
    design    = design_formula
)
mcols(dds)$gene_name <- gene_meta$gene_name

message("Running DESeq()...")
dds <- DESeq(dds, quiet = FALSE)

message("Saving fitted DDS object: ", paste0(opt$prefix, "_dds.rds"))
saveRDS(dds, file = paste0(opt$prefix, "_dds.rds"))

# ---- VST normalized matrix --------------------------------------------------
message("Computing VST (blind = FALSE)...")
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_mat <- assay(vst)

vst_df <- as.data.frame(vst_mat)
vst_df <- cbind(gene_id = rownames(vst_df),
                gene_name = gene_meta$gene_name,
                vst_df)
out_vst <- paste0(opt$prefix, "_vst_counts.tsv")
write.table(vst_df, file = out_vst, sep = "\t", quote = FALSE, row.names = FALSE)
message("VST matrix saved: ", out_vst)

# ---- size factors -----------------------------------------------------------
sf_df <- data.frame(sample = names(sizeFactors(dds)),
                    size_factor = round(sizeFactors(dds), 6))
write.csv(sf_df, file = paste0(opt$prefix, "_size_factors.csv"),
          row.names = FALSE, quote = FALSE)

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
writeLines(
    c("DESEQ2_NORMALIZATION:",
      paste0("    R: ", ver$R.version$version.string),
      paste0("    DESeq2: ", as.character(packageVersion("DESeq2"))),
      paste0("    BiocParallel: ", as.character(packageVersion("BiocParallel")))),
    "versions.yml"
)
message("Normalization complete.")
