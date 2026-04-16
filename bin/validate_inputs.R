#!/usr/bin/env Rscript
# validate_inputs.R — validate samplesheet, contrasts, and counts matrix

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option("--samplesheet", type = "character", help = "Path to samplesheet.csv"),
    make_option("--contrasts",   type = "character", help = "Path to contrasts.csv"),
    make_option("--counts",      type = "character", help = "Path to counts matrix TSV")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- helpers ----------------------------------------------------------------
stop_if <- function(cond, msg) if (cond) stop(msg, call. = FALSE)
warn_if <- function(cond, msg) if (cond) warning(msg, call. = FALSE)

# ---- load inputs ------------------------------------------------------------
message("Loading samplesheet: ", opt$samplesheet)
ss <- tryCatch(
    read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Cannot read samplesheet: ", e$message, call. = FALSE)
)

message("Loading contrasts: ", opt$contrasts)
ct <- tryCatch(
    read.csv(opt$contrasts, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Cannot read contrasts: ", e$message, call. = FALSE)
)

message("Loading counts matrix (header only): ", opt$counts)
counts_header <- tryCatch(
    colnames(read.delim(opt$counts, nrows = 1, check.names = FALSE)),
    error = function(e) stop("Cannot read counts matrix: ", e$message, call. = FALSE)
)

# ---- validate samplesheet ---------------------------------------------------
message("Validating samplesheet...")

required_ss_cols <- c("sample", "condition")
missing_ss <- setdiff(required_ss_cols, colnames(ss))
stop_if(length(missing_ss) > 0,
    paste0("Samplesheet missing required columns: ", paste(missing_ss, collapse = ", ")))

stop_if(nrow(ss) < 2, "Samplesheet must contain at least 2 samples.")

stop_if(anyDuplicated(ss$sample) > 0,
    paste0("Duplicate sample names in samplesheet: ",
           paste(ss$sample[duplicated(ss$sample)], collapse = ", ")))

stop_if(any(is.na(ss$sample) | ss$sample == ""),
    "Samplesheet contains empty or NA sample names.")
stop_if(any(is.na(ss$condition) | ss$condition == ""),
    "Samplesheet contains empty or NA condition values.")

if ("batch" %in% colnames(ss)) {
    stop_if(any(is.na(ss$batch) | ss$batch == ""),
        "Batch column contains empty or NA values. Remove the column if batch is not applicable.")
    message("  Batch column detected: design will use ~ batch + condition")
} else {
    message("  No batch column: design will use ~ condition")
}

# ---- validate contrasts -----------------------------------------------------
message("Validating contrasts...")

required_ct_cols <- c("contrast_id", "variable", "reference", "treatment")
missing_ct <- setdiff(required_ct_cols, colnames(ct))
stop_if(length(missing_ct) > 0,
    paste0("Contrasts file missing required columns: ", paste(missing_ct, collapse = ", ")))

stop_if(nrow(ct) < 1, "Contrasts file must contain at least one contrast.")

stop_if(anyDuplicated(ct$contrast_id) > 0,
    paste0("Duplicate contrast_id values: ",
           paste(ct$contrast_id[duplicated(ct$contrast_id)], collapse = ", ")))

# check that all referenced condition levels exist in samplesheet
for (i in seq_len(nrow(ct))) {
    var_col <- ct$variable[i]
    stop_if(!var_col %in% colnames(ss),
        sprintf("Contrast '%s': variable '%s' not found in samplesheet columns.",
                ct$contrast_id[i], var_col))
    available_levels <- unique(ss[[var_col]])
    stop_if(!ct$reference[i] %in% available_levels,
        sprintf("Contrast '%s': reference '%s' not found in samplesheet column '%s'.",
                ct$contrast_id[i], ct$reference[i], var_col))
    stop_if(!ct$treatment[i] %in% available_levels,
        sprintf("Contrast '%s': treatment '%s' not found in samplesheet column '%s'.",
                ct$contrast_id[i], ct$treatment[i], var_col))
    stop_if(ct$reference[i] == ct$treatment[i],
        sprintf("Contrast '%s': reference and treatment are the same ('%s').",
                ct$contrast_id[i], ct$reference[i]))
}

# ---- validate counts vs samplesheet -----------------------------------------
message("Validating counts matrix columns vs samplesheet samples...")

counts_samples <- counts_header[!counts_header %in% c("gene_id", "gene_name")]
ss_samples     <- ss$sample

missing_in_counts <- setdiff(ss_samples, counts_samples)
stop_if(length(missing_in_counts) > 0,
    paste0("Samples in samplesheet not found in counts matrix: ",
           paste(missing_in_counts, collapse = ", ")))

extra_in_counts <- setdiff(counts_samples, ss_samples)
warn_if(length(extra_in_counts) > 0,
    paste0("Counts matrix has extra columns not in samplesheet (will be ignored): ",
           paste(extra_in_counts, collapse = ", ")))

# ---- summary ----------------------------------------------------------------
message(sprintf(
    "Validation passed: %d samples, %d condition levels, %d contrasts.",
    nrow(ss),
    length(unique(ss$condition)),
    nrow(ct)
))

# ---- write validated copies -------------------------------------------------
write.csv(ss, "validated_samplesheet.csv", row.names = FALSE, quote = FALSE)
write.csv(ct, "validated_contrasts.csv",   row.names = FALSE, quote = FALSE)

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
writeLines(
    c("VALIDATE_INPUTS:",
      paste0("    R: ", ver$R.version$version.string),
      paste0("    optparse: ", as.character(packageVersion("optparse")))),
    "versions.yml"
)
message("Done.")
