#!/usr/bin/env Rscript
# regenerate_gsea_csvs.R — Re-export GSEA CSVs with correct quoting
#
# For each contrast in the GSEA output directory:
#   1. Re-export GO_all_gsea.csv from the gseaResult RDS objects (fixes leading_edge quoting bug)
#   2. Re-run MSigDB Hallmarks GSEA from the DGE CSV (fixes fgsea seed= bug)
#   3. Update the RDS with the hallmarks object
#
# Usage:
#   Rscript bin/regenerate_gsea_csvs.R --gsea_dir demo_results/gsea --dge_dir demo_results/dge --organism human

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(fgsea)
    library(msigdbr)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
})

option_list <- list(
    make_option("--gsea_dir",  type = "character", help = "Path to GSEA output directory (parent of per-contrast dirs)"),
    make_option("--dge_dir",   type = "character", help = "Path to DGE output directory (parent of per-contrast dirs)"),
    make_option("--organism",  type = "character", default = "human", help = "human or mouse [default: human]"),
    make_option("--min_gs",    type = "integer",   default = 15),
    make_option("--max_gs",    type = "integer",   default = 500),
    make_option("--n_perm",    type = "integer",   default = 1000)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$gsea_dir) || is.null(opt$dge_dir)) stop("--gsea_dir and --dge_dir are required.")

org_db <- if (opt$organism == "human") org.Hs.eg.db else org.Mm.eg.db

contrast_dirs <- list.dirs(opt$gsea_dir, full.names = TRUE, recursive = FALSE)
message(sprintf("Found %d contrast directories.", length(contrast_dirs)))

for (cdir in contrast_dirs) {
    contrast_id <- basename(cdir)
    message("\n====== Processing: ", contrast_id, " ======")

    rds_path <- file.path(cdir, paste0(contrast_id, "_gsea_dashboard_data.rds"))
    if (!file.exists(rds_path)) { message("  RDS not found, skipping."); next }

    gsea_objects <- readRDS(rds_path)

    # ---- 1. Regenerate GO CSVs from RDS objects ---------------------------
    message("  Regenerating GO CSVs from RDS...")
    go_results <- list()
    for (ont in c("BP", "MF", "CC")) {
        key <- paste0("go_", tolower(ont))
        obj <- gsea_objects[[key]]
        if (is.null(obj)) { message("    GO ", ont, ": NULL in RDS, skipping."); next }
        df <- as.data.frame(obj)
        df$ontology <- ont
        go_results[[ont]] <- df
        # Write per-ontology CSV (with quoting — write.csv default)
        out_f <- file.path(cdir, paste0(contrast_id, "_GO_", ont, "_gsea.csv"))
        write.csv(df, file = out_f, row.names = FALSE)
        message("    GO ", ont, ": wrote ", nrow(df), " rows → ", basename(out_f))
    }
    if (length(go_results) > 0) {
        go_all <- bind_rows(go_results)
        out_f  <- file.path(cdir, paste0(contrast_id, "_GO_all_gsea.csv"))
        write.csv(go_all, file = out_f, row.names = FALSE)
        message("  GO combined: wrote ", nrow(go_all), " rows → ", basename(out_f))
    }

    # ---- 2. Re-run Hallmarks if NULL in RDS --------------------------------
    if (is.null(gsea_objects$hallmarks)) {
        message("  Hallmarks is NULL — re-running fgsea...")

        # Find DGE CSV — try several naming conventions
        dge_csv <- file.path(opt$dge_dir, contrast_id, paste0(contrast_id, "_DESeq2_results.csv"))
        if (!file.exists(dge_csv)) {
            message("  DGE CSV not found at: ", dge_csv, " — trying alternative paths...")
            candidates <- Sys.glob(file.path(opt$dge_dir, contrast_id, "*results*.csv"))
            dge_csv <- if (length(candidates) > 0) candidates[1] else NA_character_
        }
        if (is.na(dge_csv) || !file.exists(dge_csv)) {
            message("  DGE CSV not found, skipping Hallmarks.")
        } else {
            message("  Loading DGE from: ", dge_csv)
            dge <- read.csv(dge_csv, stringsAsFactors = FALSE, check.names = FALSE)
            dge <- dge[!is.na(dge$pvalue) & !is.na(dge$log2FoldChange), ]
            dge$rank_metric <- sign(dge$log2FoldChange) * (-log10(dge$pvalue + .Machine$double.eps))

            # Map Ensembl → Entrez
            entrez_ids <- tryCatch(
                mapIds(org_db, keys = dge$gene_id, column = "ENTREZID",
                       keytype = "ENSEMBL", multiVals = "first"),
                error = function(e) { message("  ID mapping failed: ", e$message); NULL }
            )
            if (is.null(entrez_ids)) next

            dge$entrez_id <- entrez_ids[dge$gene_id]
            dge_mapped <- dge[!is.na(dge$entrez_id), ]
            ranked_entrez <- setNames(dge_mapped$rank_metric, dge_mapped$entrez_id)
            ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)
            ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]

            tryCatch({
                msig_species <- if (opt$organism == "human") "Homo sapiens" else "Mus musculus"
                h_sets <- msigdbr(species = msig_species, category = "H")
                h_list <- split(h_sets$entrez_gene, h_sets$gs_name)
                h_list <- lapply(h_list, function(x) as.character(unique(x)))

                set.seed(42)  # use set.seed() — fgsea() has no seed= param in v1.32+
                res_hall <- fgsea(pathways    = h_list,
                                  stats       = ranked_entrez,
                                  minSize     = opt$min_gs,
                                  maxSize     = opt$max_gs,
                                  nPermSimple = opt$n_perm,
                                  eps         = 0)
                res_hall$padj <- p.adjust(res_hall$pval, method = "BH")

                # Save native object in RDS (before collapsing leadingEdge)
                gsea_objects$hallmarks <- data.table::copy(res_hall)

                # Collapse leadingEdge list column → character for CSV
                res_hall_csv <- as.data.frame(res_hall)
                res_hall_csv$leadingEdge <- sapply(res_hall$leadingEdge, paste, collapse = "/")

                out_f <- file.path(cdir, paste0(contrast_id, "_Hallmarks_gsea.csv"))
                write.csv(res_hall_csv, file = out_f, row.names = FALSE)
                message("  Hallmarks: wrote ", nrow(res_hall_csv), " rows → ", basename(out_f))
            }, error = function(e) {
                message("  Hallmarks GSEA failed: ", e$message)
            })
        }
    } else {
        message("  Hallmarks already present in RDS, skipping re-run.")
    }

    # ---- 3. Save updated RDS -----------------------------------------------
    saveRDS(gsea_objects, file = rds_path)
    message("  Updated RDS saved: ", basename(rds_path))
}

message("\nDone. All contrasts processed.")
