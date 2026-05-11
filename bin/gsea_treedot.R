#!/usr/bin/env Rscript
# gsea_treedot.R — Cross-contrast GSEA comparison: TreeDot plot
# Produces a dendrogram + dot-heatmap comparing GSEA enrichment across all contrasts.
# Input: all *_DESeq2_results.csv files staged in the Nextflow work directory.

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(clusterProfiler)
    library(ReactomePA)
    library(enrichplot)
    library(AnnotationDbi)
    library(ggtree)
    library(cowplot)
    library(RColorBrewer)
})

option_list <- list(
    make_option("--organism",   type = "character", default = "human",
                help = "Organism: human or mouse [default: %default]"),
    make_option("--databases",  type = "character", default = "GO:BP,KEGG,Reactome",
                help = "Comma-separated databases: GO:BP, GO:MF, GO:CC, KEGG, Reactome [default: %default]"),
    make_option("--top_paths",  type = "integer",   default = 5L,
                help = "Top N pathways per contrast shown in the plot [default: %default]"),
    make_option("--clust_num",  type = "integer",   default = 3L,
                help = "Number of dendrogram clusters [default: %default]"),
    make_option("--ora_type",   type = "character", default = "GO",
                help = "ORA database for branch labels: GO or KEGG [default: %default]"),
    make_option("--ora_ont",    type = "character", default = "BP",
                help = "ORA GO ontology [default: %default]"),
    make_option("--ora_min_gs", type = "integer",   default = 10L,
                help = "ORA minimum gene set size [default: %default]"),
    make_option("--ora_max_gs", type = "integer",   default = 500L,
                help = "ORA maximum gene set size [default: %default]"),
    make_option("--ora_padj",   type = "double",    default = 1.0,
                help = "ORA adj.p cutoff (1 = include all terms) [default: %default]"),
    make_option("--min_gs",     type = "integer",   default = 15L,
                help = "GSEA minimum gene set size [default: %default]"),
    make_option("--max_gs",     type = "integer",   default = 500L,
                help = "GSEA maximum gene set size [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Organism setup ---------------------------------------------------------
if (opt$organism == "human") {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    org_db     <- org.Hs.eg.db
    kegg_org   <- "hsa"
    org_db_str <- "org.Hs.eg.db"
} else {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    org_db     <- org.Mm.eg.db
    kegg_org   <- "mmu"
    org_db_str <- "org.Mm.eg.db"
}

# ---- Discover DGE CSV files staged in the work directory --------------------
dge_files <- sort(Sys.glob("*_DESeq2_results.csv"))

if (length(dge_files) < 2) {
    msg <- if (length(dge_files) == 0) "No *_DESeq2_results.csv files found."
           else "Only 1 contrast found — TreeDot requires >= 2 contrasts."
    message(msg, " Skipping TreeDot.")
    writeLines(c("GSEA_TREEDOT:",
                 paste0("    status: skipped (", length(dge_files), " contrast(s))")),
               "versions.yml")
    quit(save = "no", status = 0)
}

contrast_ids <- sub("_DESeq2_results\\.csv$", "", basename(dge_files))
message("Building TreeDot for ", length(dge_files), " contrasts:")
for (cid in contrast_ids) message("  ", cid)

# ---- Build ranked gene lists ------------------------------------------------
message("Building ranked gene lists (metric: sign(log2FC) x -log10(pvalue))...")
ranked_ensembl_list <- list()
ranked_entrez_list  <- list()

for (i in seq_along(dge_files)) {
    cid <- contrast_ids[i]
    df  <- read.csv(dge_files[i], stringsAsFactors = FALSE, check.names = FALSE)
    df  <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), ]
    df$rank_metric <- sign(df$log2FoldChange) * (-log10(df$pvalue + .Machine$double.eps))

    # ENSEMBL ranked list (for GO with ENSEMBL keys)
    rl_ens <- setNames(df$rank_metric, df$gene_id)
    rl_ens <- sort(rl_ens, decreasing = TRUE)
    ranked_ensembl_list[[cid]] <- rl_ens[!duplicated(names(rl_ens))]

    # ENTREZ ranked list (for KEGG)
    ez <- mapIds(org_db, keys = df$gene_id, column = "ENTREZID",
                 keytype = "ENSEMBL", multiVals = "first")
    df$entrez_id <- ez[df$gene_id]
    df_ez <- df[!is.na(df$entrez_id), ]
    if (nrow(df_ez) > 0) {
        rl_ez <- setNames(df_ez$rank_metric, df_ez$entrez_id)
        rl_ez <- sort(rl_ez, decreasing = TRUE)
        ranked_entrez_list[[cid]] <- rl_ez[!duplicated(names(rl_ez))]
    }
}

# ---- compareCluster ---------------------------------------------------------
message("Running compareCluster (database=", opt$database, ", go_ont=", opt$go_ont, ")...")
set.seed(42)

cmp <- tryCatch({
    if (opt$database == "GO") {
        compareCluster(
            geneClusters = ranked_ensembl_list,
            fun          = "gseGO",
            OrgDb        = org_db,
            keyType      = "ENSEMBL",
            ont          = opt$go_ont,
            pvalueCutoff = 1,
            minGSSize    = opt$min_gs,
            maxGSSize    = opt$max_gs,
            eps          = 0,
            seed         = TRUE,
            verbose      = FALSE
        )
    } else {
        valid_el <- ranked_entrez_list[lengths(ranked_entrez_list) > 0]
        if (length(valid_el) < 2) stop("< 2 contrasts have Entrez IDs — cannot run KEGG TreeDot.")
        compareCluster(
            geneClusters = valid_el,
            fun          = "gseKEGG",
            organism     = kegg_org,
            pvalueCutoff = 1,
            minGSSize    = opt$min_gs,
            maxGSSize    = opt$max_gs,
            eps          = 0,
            seed         = TRUE,
            verbose      = FALSE
        )
    }
}, error = function(e) { message("compareCluster error: ", e$message); NULL })

if (is.null(cmp) || nrow(as.data.frame(cmp)) == 0) {
    message("No enriched terms found — skipping TreeDot plot.")
    writeLines(c("GSEA_TREEDOT:", "    status: skipped (no enriched terms)"), "versions.yml")
    quit(save = "no", status = 0)
}
message("compareCluster: ", nrow(as.data.frame(cmp)), " term-contrast associations found.")

# ---- pairwise_termsim -------------------------------------------------------
# Required for the dendrogram: fills the @termsim slot of the compareClusterResult.
message("Computing pairwise term similarity...")
cmp <- tryCatch(
    pairwise_termsim(cmp),
    error = function(e) { message("pairwise_termsim warning: ", e$message); cmp }
)

# Save RDS — dashboard reuses this to replot without re-running compareCluster
saveRDS(cmp, "gsea_treedot.rds")
message("Saved: gsea_treedot.rds")

# ---- TreeDot plot function --------------------------------------------------
# Adapted from scPlot_treedot for bulk contrasts.
# cmp      : compareClusterResult with pairwise_termsim() applied
# keytype  : gene ID type in core_enrichment slot ("ENSEMBL" for GO, "ENTREZID" for KEGG)
plot_treedot <- function(cmp, top_paths = 5, clust_num = 3,
                         ora_type = "GO", ora_ont = "BP",
                         ora_min_gs = 10, ora_max_gs = 500, ora_padj = 1,
                         org_db_str = "org.Hs.eg.db", kegg_org = "hsa",
                         keytype = "ENSEMBL") {

    # Fortify — top N pathways per contrast
    fort <- fortify(cmp, showCategory = top_paths, includeAll = TRUE, split = NULL)
    fort$Cluster <- sub("\n.*", "", fort$Cluster)
    fort$geneID  <- fort$core_enrichment
    if (nrow(fort) == 0) stop("No pathways after fortify — increase top_paths or check GSEA results.")

    # Build pathway × contrast presence matrix (used for termsim indexing)
    id_unique <- unique(fort$Description)
    cl_unique <- unique(fort$Cluster)
    mat <- matrix(0, nrow = length(id_unique), ncol = length(cl_unique),
                  dimnames = list(id_unique, cl_unique))
    for (i in seq_len(nrow(fort))) mat[fort$Description[i], fort$Cluster[i]] <- 1
    id_mat <- as.data.frame(mat)

    # Term similarity → hierarchical clustering
    fill_termsim <- function(x, keep) {
        ts <- x@termsim[keep, keep, drop = FALSE]
        ts[is.na(ts)] <- 0
        ts2 <- ts + t(ts)
        diag(ts2) <- 1
        ts2
    }
    ts2 <- tryCatch(
        fill_termsim(cmp, rownames(id_mat)),
        error = function(e) {
            message("Termsim unavailable, using identity matrix: ", e$message)
            diag_m <- diag(nrow(id_mat))
            rownames(diag_m) <- colnames(diag_m) <- rownames(id_mat)
            diag_m
        }
    )
    hc        <- stats::hclust(stats::as.dist(1 - ts2), method = "ward.D")
    clust_num <- min(clust_num, nrow(id_mat))
    clus      <- stats::cutree(hc, clust_num)

    # ORA per cluster → branch label (top term from ORA as cluster name)
    keywords <- character(clust_num)
    for (i in seq_len(clust_num)) {
        paths_i    <- names(clus[clus == i])
        cluster_df <- fort[fort$Description %in% paths_i, ]
        genes_i    <- unique(unlist(str_split(cluster_df$geneID, "/")))
        genes_i    <- genes_i[nzchar(genes_i) & !is.na(genes_i)]

        # enrichKEGG requires ENTREZ IDs; convert if the main GSEA used ENSEMBL keys
        genes_for_ora <- genes_i
        if (ora_type == "KEGG" && keytype != "ENTREZID") {
            genes_for_ora <- tryCatch(
                na.omit(as.character(mapIds(get(org_db_str), keys = genes_i,
                                           column = "ENTREZID", keytype = keytype,
                                           multiVals = "first"))),
                error = function(e) { message("ID conversion for KEGG ORA: ", e$message); character(0) }
            )
        }

        ora_res <- tryCatch({
            if (ora_type == "GO") {
                enrichGO(gene = genes_for_ora, OrgDb = org_db_str, keyType = keytype,
                         ont = ora_ont, pvalueCutoff = ora_padj, pAdjustMethod = "BH",
                         qvalueCutoff = 1, minGSSize = ora_min_gs, maxGSSize = ora_max_gs)
            } else {
                enrichKEGG(gene = genes_for_ora, organism = kegg_org,
                           pvalueCutoff = ora_padj, pAdjustMethod = "BH",
                           qvalueCutoff = 1, minGSSize = ora_min_gs, maxGSSize = ora_max_gs)
            }
        }, error = function(e) { message("ORA cluster ", i, ": ", e$message); NULL })

        keywords[i] <- if (!is.null(ora_res) && nrow(ora_res@result) > 0 &&
                            !is.na(ora_res@result$Description[1]))
            ora_res@result$Description[1]
        else
            paste0("Cluster ", i)
    }

    kw_maxlen <- max(nchar(keywords), 17)

    # ggtree dendrogram coloured by clade / cluster
    g_split     <- split(names(clus), clus)
    p_tree      <- ggtree(hc, size = 1.2)
    clades      <- sapply(g_split, function(n) MRCA(p_tree, n))
    p_tree      <- groupClade(p_tree, clades, group_name = "SubTree_ORA") +
        aes(color = SubTree_ORA)
    pal         <- c(brewer.pal(8, "Dark2"), brewer.pal(6, "Set1"))

    ggtree_full  <- p_tree +
        scale_color_manual(name = "SubTree ORA",
                           values = pal, breaks = seq_len(clust_num), labels = keywords) +
        theme(legend.position = "right", legend.justification = c(0, 1.5))
    ggtree_noleg <- ggtree_full + theme(legend.position = "none")

    # Prepare dot plot data
    fort$log_p.adjust <- -log10(fort$p.adjust)
    fort$Description  <- factor(fort$Description, levels = hc$labels[hc$order])
    cl_levels <- if (is.factor(cmp@compareClusterResult$Cluster))
        sub("\n.*", "", levels(cmp@compareClusterResult$Cluster))
    else
        sub("\n.*", "", unique(as.character(cmp@compareClusterResult$Cluster)))
    fort$Cluster <- factor(fort$Cluster, levels = cl_levels)

    dotplot_full <- fort %>%
        ggplot(aes(x = Cluster, y = Description, color = NES, size = log_p.adjust)) +
        geom_point() +
        scale_y_discrete(position = "right") +
        scale_color_gradient2(low = "blue4", mid = "white", high = "red") +
        theme_cowplot(font_size = 14) +
        theme(axis.line   = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13),
              axis.text.y = element_text(size = 12)) +
        xlab("Contrast") + ylab("") +
        guides(size  = guide_legend(title = "-log10(p.adj)"),
               color = guide_colorbar(title = "NES")) +
        theme(axis.ticks    = element_blank(),
              legend.position  = "right", legend.justification = c(0, 0),
              legend.text      = element_text(size = 12),
              legend.title     = element_text(size = 13))
    dotplot_noleg <- dotplot_full + theme(legend.position = "none")

    # Assemble final plot: dendrogram | dot-heatmap | legends
    leg_tree <- get_legend(ggtree_full)
    leg_dot  <- get_legend(dotplot_full)
    row_main <- plot_grid(ggtree_noleg, NULL, dotplot_noleg,
                          nrow = 1, rel_widths = c(0.3, -0.05, 2), align = "h")
    leg_col  <- plot_grid(leg_dot, NULL, leg_tree,
                          ncol = 1, rel_heights = c(1, -0.5, 1))
    plot_grid(row_main, leg_col, NULL,
              nrow = 1, rel_widths = c(1, 0.1, kw_maxlen * 0.006))
}

# ---- versions ---------------------------------------------------------------
ver      <- sessionInfo()
pkg_list <- c("clusterProfiler", "ReactomePA", "enrichplot", "ggtree", "cowplot",
              "RColorBrewer", "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi",
              "ggplot2", "dplyr")
ver_lines <- c("GSEA_TREEDOT:", paste0("    R: ", ver$R.version$version.string))
for (pkg in pkg_list) {
    tryCatch(ver_lines <<- c(ver_lines, paste0("    ", pkg, ": ", packageVersion(pkg))),
             error = function(e) NULL)
}
writeLines(ver_lines, "versions.yml")
message("TreeDot complete.")
