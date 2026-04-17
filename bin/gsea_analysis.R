#!/usr/bin/env Rscript
# gsea_analysis.R — Gene set enrichment analysis for a single contrast
# Runs GO (BP/MF/CC), KEGG, MSigDB Hallmarks, Reactome with clusterProfiler/fgsea.

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
    library(clusterProfiler)
    library(fgsea)
    library(msigdbr)
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
})

# Container for native result objects (gseaResult / data.table) — saved as RDS
# for the dashboard to use enrichplot functions directly.
gsea_objects <- list(
    go_bp    = NULL,
    go_mf    = NULL,
    go_cc    = NULL,
    kegg     = NULL,
    hallmarks = NULL,
    reactome = NULL
)

option_list <- list(
    make_option("--dge_results", type = "character", help = "DESeq2 results CSV (gene_id, log2FoldChange, pvalue, padj)"),
    make_option("--contrast_id", type = "character", help = "Contrast identifier"),
    make_option("--organism",    type = "character", default = "human",
                help = "Organism: human or mouse [default: human]"),
    make_option("--min_gs",      type = "integer",   default = 15,
                help = "Min gene set size [default: 15]"),
    make_option("--max_gs",      type = "integer",   default = 500,
                help = "Max gene set size [default: 500]"),
    make_option("--pval_cutoff", type = "double",    default = 0.05,
                help = "Adjusted p-value cutoff for plots [default: 0.05]"),
    make_option("--n_perm",      type = "integer",   default = 1000,
                help = "Number of fgsea permutations [default: 1000]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- setup ------------------------------------------------------------------
org_db  <- if (opt$organism == "human") org.Hs.eg.db else org.Mm.eg.db
kegg_org <- if (opt$organism == "human") "hsa" else "mmu"
react_org <- if (opt$organism == "human") "human" else "mouse"
contrast_id <- opt$contrast_id

message("Running GSEA for contrast: ", contrast_id)
message("Organism: ", opt$organism)

# ---- load DGE results -------------------------------------------------------
message("Loading DGE results: ", opt$dge_results)
dge <- read.csv(opt$dge_results, stringsAsFactors = FALSE, check.names = FALSE)

# remove rows with NA pvalue or log2FC
dge <- dge[!is.na(dge$pvalue) & !is.na(dge$log2FoldChange), ]

# ---- ranking metric: sign(log2FC) * -log10(pvalue) --------------------------
dge$rank_metric <- sign(dge$log2FoldChange) * (-log10(dge$pvalue + .Machine$double.eps))

# ---- Ensembl → Entrez ID mapping -------------------------------------------
message("Mapping Ensembl → Entrez IDs...")
ensembl_ids <- dge$gene_id
entrez_ids  <- mapIds(org_db,
                      keys      = ensembl_ids,
                      column    = "ENTREZID",
                      keytype   = "ENSEMBL",
                      multiVals = "first")

n_mapped <- sum(!is.na(entrez_ids))
message(sprintf("  Mapped %d / %d genes to Entrez IDs.", n_mapped, length(ensembl_ids)))
if (n_mapped < 100) warning("Fewer than 100 genes mapped — GSEA results may be unreliable.", call. = FALSE)

dge$entrez_id <- entrez_ids[dge$gene_id]
dge_mapped <- dge[!is.na(dge$entrez_id), ]

# ranked gene list (Entrez)
ranked_entrez <- setNames(dge_mapped$rank_metric, dge_mapped$entrez_id)
ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)
ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]

# ranked gene list (Ensembl — for GO with ENSEMBL keys)
ranked_ensembl <- setNames(dge$rank_metric, dge$gene_id)
ranked_ensembl <- sort(ranked_ensembl, decreasing = TRUE)
ranked_ensembl <- ranked_ensembl[!duplicated(names(ranked_ensembl))]

# ---- helper: save ggplot as PDF + PNG --------------------------------------
save_plot <- function(base, p, width = 10, height = 7) {
    for (ext in c("pdf", "png")) {
        f <- paste0(base, ".", ext)
        if (ext == "png") png(f, width = width * 100, height = height * 100, res = 100)
        else              pdf(f, width = width, height = height)
        print(p)
        dev.off()
    }
}

# ---- helper: dot plot for enrichment results --------------------------------
make_dotplot <- function(df, title, top_n = 20) {
    # clusterProfiler uses 'p.adjust'; fgsea uses 'padj' — normalise
    if (!"padj" %in% colnames(df) && "p.adjust" %in% colnames(df)) {
        df$padj <- df$p.adjust
    }
    df <- df %>%
        filter(padj < opt$pval_cutoff) %>%
        arrange(padj) %>%
        head(top_n) %>%
        mutate(Description = stringr::str_wrap(Description, width = 50),
               sign        = ifelse(NES > 0, "Up", "Down"))
    if (nrow(df) == 0) return(NULL)
    ggplot(df, aes(x = NES, y = reorder(Description, NES),
                   colour = padj, size = size)) +
        geom_point() +
        scale_colour_gradient(low = "#d73027", high = "#4575b4",
                              name = "adj. p-value") +
        scale_size_continuous(name = "Gene set size", range = c(3, 10)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
        labs(title = title, x = "Normalised Enrichment Score", y = NULL) +
        theme_bw(base_size = 11) +
        theme(axis.text.y = element_text(size = 9))
}

# ---- GO enrichment (BP, MF, CC) ---------------------------------------------
message("Running GO GSEA (BP, MF, CC)...")
go_results <- list()
for (ont in c("BP", "MF", "CC")) {
    message("  GO: ", ont)
    set.seed(42)
    res_go <- tryCatch(
        gseGO(geneList     = ranked_ensembl,
              OrgDb        = org_db,
              ont          = ont,
              keyType      = "ENSEMBL",
              minGSSize    = opt$min_gs,
              maxGSSize    = opt$max_gs,
              pvalueCutoff = 1,
              verbose      = FALSE,
              eps          = 0,
              seed         = TRUE),
        error = function(e) { message("    GO ", ont, " failed: ", e$message); NULL }
    )
    if (!is.null(res_go)) {
        gsea_objects[[paste0("go_", tolower(ont))]] <- res_go  # save gseaResult object (for loop: no new env, use <-)
        df <- as.data.frame(res_go)
        df$ontology <- ont
        go_results[[ont]] <- df
        out_f <- paste0(contrast_id, "_GO_", ont, "_gsea.csv")
        write.csv(df, file = out_f, row.names = FALSE)

        df_plot <- df %>%
            mutate(size = sapply(strsplit(core_enrichment, "/"), length))
        p <- make_dotplot(df_plot, paste0("GO ", ont, " — ", contrast_id))
        if (!is.null(p)) save_plot(paste0(contrast_id, "_GO_", ont, "_dotplot"), p)
    }
}

# combine all GO
if (length(go_results) > 0) {
    go_all <- bind_rows(go_results)
    write.csv(go_all, file = paste0(contrast_id, "_GO_all_gsea.csv"),
              row.names = FALSE)
}

# ---- KEGG -------------------------------------------------------------------
message("Running KEGG GSEA (requires internet)...")
tryCatch({
    set.seed(42)
    res_kegg <- gseKEGG(geneList     = ranked_entrez,
                        organism     = kegg_org,
                        minGSSize    = opt$min_gs,
                        maxGSSize    = opt$max_gs,
                        pvalueCutoff = 1,
                        verbose      = FALSE,
                        eps          = 0,
                        seed         = TRUE,
                        use_internal_data = FALSE)
    gsea_objects$kegg <<- res_kegg  # save gseaResult object
    df_kegg <- as.data.frame(res_kegg)
    write.csv(df_kegg, file = paste0(contrast_id, "_KEGG_gsea.csv"),
              row.names = FALSE)
    df_kegg$size <- sapply(strsplit(df_kegg$core_enrichment, "/"), length)
    p <- make_dotplot(df_kegg, paste0("KEGG — ", contrast_id))
    if (!is.null(p)) save_plot(paste0(contrast_id, "_KEGG_dotplot"), p)
    message("  KEGG done.")
}, error = function(e) {
    message("  KEGG GSEA failed (offline or API issue): ", e$message)
    # write empty file so process doesn't fail
    write.csv(data.frame(), file = paste0(contrast_id, "_KEGG_gsea.csv"),
              row.names = FALSE)
})

# ---- MSigDB Hallmarks -------------------------------------------------------
message("Running MSigDB Hallmarks GSEA...")
tryCatch({
    msig_species <- if (opt$organism == "human") "Homo sapiens" else "Mus musculus"
    h_sets   <- msigdbr(species = msig_species, category = "H")
    h_list   <- split(h_sets$entrez_gene, h_sets$gs_name)
    h_list   <- lapply(h_list, function(x) as.character(unique(x)))

    set.seed(42)  # fgsea v1.32+ has no seed= param; use set.seed() instead
    res_hall <- fgsea(pathways    = h_list,
                      stats       = ranked_entrez,
                      minSize     = opt$min_gs,
                      maxSize     = opt$max_gs,
                      nPermSimple = opt$n_perm,
                      eps         = 0)
    res_hall$padj <- p.adjust(res_hall$pval, method = "BH")
    gsea_objects$hallmarks <<- data.table::copy(res_hall)  # save native data.table (leadingEdge as list)
    # convert to data.frame to avoid data.table list-column issues in write.csv
    res_hall_df <- as.data.frame(res_hall)
    res_hall_df$leadingEdge <- sapply(res_hall_df$leadingEdge, paste, collapse = "/")

    write.csv(res_hall_df, file = paste0(contrast_id, "_Hallmarks_gsea.csv"),
              row.names = FALSE)

    hall_plot <- res_hall %>%
        rename(Description = pathway, NES = NES, size = size) %>%
        mutate(size = as.numeric(size))
    p <- make_dotplot(hall_plot, paste0("MSigDB Hallmarks — ", contrast_id))
    if (!is.null(p)) save_plot(paste0(contrast_id, "_Hallmarks_dotplot"), p)
    message("  Hallmarks done.")
}, error = function(e) {
    message("  Hallmarks GSEA failed: ", e$message)
    write.csv(data.frame(), file = paste0(contrast_id, "_Hallmarks_gsea.csv"), row.names = FALSE)
})

# ---- Reactome ---------------------------------------------------------------
message("Running Reactome GSEA (requires internet)...")
tryCatch({
    set.seed(42)
    res_react <- gsePathway(geneList     = ranked_entrez,
                            organism     = react_org,
                            minGSSize    = opt$min_gs,
                            maxGSSize    = opt$max_gs,
                            pvalueCutoff = 1,
                            verbose      = FALSE,
                            eps          = 0,
                            seed         = TRUE)
    gsea_objects$reactome <<- res_react  # save gseaResult object
    df_react <- as.data.frame(res_react)
    write.csv(df_react, file = paste0(contrast_id, "_Reactome_gsea.csv"),
              row.names = FALSE)
    df_react$size <- sapply(strsplit(df_react$core_enrichment, "/"), length)
    p <- make_dotplot(df_react, paste0("Reactome — ", contrast_id))
    if (!is.null(p)) save_plot(paste0(contrast_id, "_Reactome_dotplot"), p)
    message("  Reactome done.")
}, error = function(e) {
    message("  Reactome GSEA failed (offline or API issue): ", e$message)
    write.csv(data.frame(), file = paste0(contrast_id, "_Reactome_gsea.csv"), row.names = FALSE)
})

# ---- save dashboard RDS (native enrichResult/gseResult objects) -------------
rds_path <- paste0(contrast_id, "_gsea_dashboard_data.rds")
saveRDS(gsea_objects, file = rds_path)
message("Saved dashboard RDS: ", rds_path)

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
pkg_versions <- c("clusterProfiler", "fgsea", "msigdbr", "ReactomePA",
                  "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi",
                  "ggplot2", "dplyr")
ver_lines <- c("GSEA_ANALYSIS:",
               paste0("    R: ", ver$R.version$version.string))
for (pkg in pkg_versions) {
    tryCatch(
        ver_lines <<- c(ver_lines, paste0("    ", pkg, ": ", as.character(packageVersion(pkg)))),
        error = function(e) NULL
    )
}
writeLines(ver_lines, "versions.yml")
message("GSEA complete: ", contrast_id)
