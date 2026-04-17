#!/usr/bin/env Rscript
# gsea_analysis.R — Gene set enrichment analysis for a single contrast
# Runs GO (BP/MF/CC), KEGG, Reactome with clusterProfiler/ReactomePA.

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
    library(ggridges)
    library(clusterProfiler)
    library(enrichplot)
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
})

# Container for native result objects (gseaResult S4) — saved as RDS
# for the dashboard to use enrichplot functions directly.
gsea_objects <- list(
    go_bp    = NULL,
    go_mf    = NULL,
    go_cc    = NULL,
    kegg     = NULL,
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
                help = "Adjusted p-value cutoff for plots [default: 0.05]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- setup ------------------------------------------------------------------
org_db    <- if (opt$organism == "human") org.Hs.eg.db else org.Mm.eg.db
kegg_org  <- if (opt$organism == "human") "hsa" else "mmu"
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
    # clusterProfiler uses 'p.adjust'; normalise to 'padj'
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

# ---- helper: ridge plot for gseaResult S4 objects (GO / KEGG / Reactome) ----
make_ridgeplot_s4 <- function(gse_obj, title, top_n = 20) {
    df <- tryCatch(as.data.frame(gse_obj), error = function(e) data.frame())
    if (nrow(df) == 0) return(NULL)
    col_p <- if ("p.adjust" %in% colnames(df)) "p.adjust" else if ("padj" %in% colnames(df)) "padj" else NULL
    n_sig <- if (!is.null(col_p)) sum(!is.na(df[[col_p]]) & df[[col_p]] < opt$pval_cutoff) else nrow(df)
    if (n_sig == 0) return(NULL)
    n_show <- min(top_n, n_sig)
    tryCatch(
        enrichplot::ridgeplot(gse_obj, showCategory = n_show) +
            labs(title = title) +
            theme_bw(base_size = 11),
        error = function(e) { message("    ridgeplot error: ", e$message); NULL }
    )
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
        gsea_objects[[paste0("go_", tolower(ont))]] <- res_go
        df <- as.data.frame(res_go)
        df$ontology <- ont
        go_results[[ont]] <- df
        write.csv(df, file = paste0(contrast_id, "_GO_", ont, "_gsea.csv"), row.names = FALSE)

        df_plot <- df %>%
            mutate(size = sapply(strsplit(core_enrichment, "/"), length))
        p <- make_dotplot(df_plot, paste0("GO ", ont, " — ", contrast_id))
        if (!is.null(p)) save_plot(paste0(contrast_id, "_GO_", ont, "_dotplot"), p)

        p_ridge <- make_ridgeplot_s4(res_go, paste0("GO ", ont, " Ridgeplot — ", contrast_id))
        if (!is.null(p_ridge)) {
            n_sig <- sum(!is.na(df$p.adjust) & df$p.adjust < opt$pval_cutoff)
            rh <- max(5, min(n_sig * 0.45 + 2, 16))
            save_plot(paste0(contrast_id, "_GO_", ont, "_ridgeplot"), p_ridge, height = rh)
        }
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
    gsea_objects$kegg <- res_kegg
    df_kegg <- as.data.frame(res_kegg)
    write.csv(df_kegg, file = paste0(contrast_id, "_KEGG_gsea.csv"),
              row.names = FALSE)
    df_kegg$size <- sapply(strsplit(df_kegg$core_enrichment, "/"), length)
    p <- make_dotplot(df_kegg, paste0("KEGG — ", contrast_id))
    if (!is.null(p)) save_plot(paste0(contrast_id, "_KEGG_dotplot"), p)

    p_ridge <- make_ridgeplot_s4(res_kegg, paste0("KEGG Ridgeplot — ", contrast_id))
    if (!is.null(p_ridge)) {
        n_sig <- sum(!is.na(df_kegg$p.adjust) & df_kegg$p.adjust < opt$pval_cutoff)
        rh <- max(5, min(n_sig * 0.45 + 2, 16))
        save_plot(paste0(contrast_id, "_KEGG_ridgeplot"), p_ridge, height = rh)
    }

    # ---- Pathview: top 10 significant KEGG pathways -------------------------
    message("  Generating Pathview for top KEGG pathways...")
    tryCatch({
        fc_entrez   <- setNames(dge_mapped$log2FoldChange, as.character(dge_mapped$entrez_id))
        df_kegg_sig <- df_kegg[order(df_kegg$p.adjust), ]
        top10_ids   <- head(df_kegg_sig$ID, 10)
        for (pw_full in top10_ids) {
            pw_num <- sub("^[a-z]+", "", pw_full)  # "hsa04110" → "04110"
            tryCatch({
                suppressMessages(
                    pathview::pathview(
                        gene.data   = fc_entrez,
                        pathway.id  = pw_num,
                        species     = kegg_org,
                        out.suffix  = contrast_id,
                        kegg.dir    = ".",
                        gene.idtype = "entrez",
                        low  = list(gene = "#4575b4"),
                        mid  = list(gene = "#ffffbf"),
                        high = list(gene = "#d73027"),
                        na.col = "grey80"
                    )
                )
                old_f <- paste0(kegg_org, pw_num, ".", contrast_id, ".png")
                new_f <- paste0(contrast_id, "_pathview_", kegg_org, pw_num, ".png")
                if (file.exists(old_f)) {
                    file.rename(old_f, new_f)
                    message("    Pathview saved: ", new_f)
                }
            }, error = function(e) message("    Pathview failed (", pw_full, "): ", e$message))
        }
    }, error = function(e) message("  Pathview block failed: ", e$message))

    message("  KEGG done.")
}, error = function(e) {
    message("  KEGG GSEA failed (offline or API issue): ", e$message)
    write.csv(data.frame(), file = paste0(contrast_id, "_KEGG_gsea.csv"),
              row.names = FALSE)
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
    gsea_objects$reactome <- res_react
    df_react <- as.data.frame(res_react)
    write.csv(df_react, file = paste0(contrast_id, "_Reactome_gsea.csv"),
              row.names = FALSE)
    df_react$size <- sapply(strsplit(df_react$core_enrichment, "/"), length)
    p <- make_dotplot(df_react, paste0("Reactome — ", contrast_id))
    if (!is.null(p)) save_plot(paste0(contrast_id, "_Reactome_dotplot"), p)

    p_ridge <- make_ridgeplot_s4(res_react, paste0("Reactome Ridgeplot — ", contrast_id))
    if (!is.null(p_ridge)) {
        n_sig <- sum(!is.na(df_react$p.adjust) & df_react$p.adjust < opt$pval_cutoff)
        rh <- max(5, min(n_sig * 0.45 + 2, 16))
        save_plot(paste0(contrast_id, "_Reactome_ridgeplot"), p_ridge, height = rh)
    }
    message("  Reactome done.")
}, error = function(e) {
    message("  Reactome GSEA failed (offline or API issue): ", e$message)
    write.csv(data.frame(), file = paste0(contrast_id, "_Reactome_gsea.csv"), row.names = FALSE)
})

# ---- save dashboard RDS (native gseaResult S4 objects) ----------------------
rds_path <- paste0(contrast_id, "_gsea_dashboard_data.rds")
saveRDS(gsea_objects, file = rds_path)
message("Saved dashboard RDS: ", rds_path)

# ---- versions ---------------------------------------------------------------
ver <- sessionInfo()
pkg_versions <- c("clusterProfiler", "enrichplot", "ReactomePA", "pathview",
                  "ggridges", "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi",
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
