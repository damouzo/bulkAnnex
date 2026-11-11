# modules/mod_dge.R — DGE tab

mod_dge_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12, h3("DGE"))
        ),
        fluidRow(
            column(3,
                card(
                    card_header("Settings"),
                    selectInput(ns("contrast"),  "Contrast", choices = NULL),
                    numericInput(ns("padj_cut"), "FDR cutoff",
                                 value = 0.05, min = 0.001, max = 0.2, step = 0.001),
                    numericInput(ns("lfc_cut"),  "|log2FC| cutoff",
                                 value = 1, min = 0, max = 10, step = 0.1),
                    numericInput(ns("n_label"),  "Label top N genes", value = 30, min = 0, max = 100),
                    sliderInput(ns("label_size"), "Label font size",
                                min = 6, max = 16, value = 9, step = 1),
                    textAreaInput(ns("must_genes"), "Genes to always label",
                                  value = "", rows = 2,
                                  placeholder = "e.g. EGR2, LTF"),
                    fluidRow(
                        column(6, checkboxInput(ns("highlight_labels"), "Highlight label", value = FALSE)),
                        column(6, checkboxInput(ns("hide_top_n"), "Hide top N", value = FALSE))
                    ),
                    actionButton(ns("apply"), "Apply settings",
                                 icon = icon("play"),
                                 class = "btn-primary w-100 mt-2"),
                    downloadButton(ns("download_volcano"), "Download PNG",
                                   class = "btn-outline-secondary w-100 mt-2"),
                    downloadButton(ns("download_volcano_pdf"), "Download PDF",
                                   class = "btn-outline-secondary w-100 mt-2")
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Volcano",
                        imageOutput(ns("volcano_plot"), width = "100%", height = "700px")
                    ),
                    nav_panel("Heatmap",
                        fluidRow(
                            column(3,
                                card(
                                    card_header("Heatmap settings"),
                                    numericInput(ns("hm_top_n"), "Top N genes (rows)",
                                                 value = 50, min = 2, max = 1000, step = 1),
                                    numericInput(ns("hm_lfc"), "|log2FC| cutoff",
                                                 value = 1, min = 0, max = 10, step = 0.1),
                                    numericInput(ns("hm_padj"), "FDR cutoff",
                                                 value = 0.05, min = 0.001, max = 0.2, step = 0.001),
                                    radioButtons(ns("hm_rank"), "Rank genes by",
                                        choices  = c("|log2FC|" = "logfc", "padj" = "padj"),
                                        selected = "logfc", inline = TRUE),
                                    radioButtons(ns("hm_mode"), "Ranking mode",
                                        choices  = c("Top N by metric" = "abs",
                                                     "Balanced (N/2 up + N/2 down)" = "balanced"),
                                        selected = "abs"),
                                    checkboxInput(ns("hm_names"), "Show gene names", value = TRUE),
                                    checkboxInput(ns("hm_cluster_rows"), "Cluster rows (genes)", value = TRUE),
                                    checkboxInput(ns("hm_cluster_cols"), "Cluster columns (samples)", value = TRUE),
                                    checkboxInput(ns("hm_dend"), "Show dendrograms", value = TRUE),
                                    checkboxInput(ns("hm_supervise"),
                                                  "Supervise columns (group by condition)", value = FALSE),
                                    actionButton(ns("hm_apply"), "Apply settings",
                                                 icon = icon("play"),
                                                 class = "btn-primary w-100 mt-2"),
                                    downloadButton(ns("hm_download"), "Download PNG",
                                                   class = "btn-outline-secondary w-100 mt-2"),
                                    downloadButton(ns("hm_download_pdf"), "Download PDF",
                                                   class = "btn-outline-secondary w-100 mt-2")
                                )
                            ),
                            column(9,
                                card(
                                    card_header("Top genes heatmap"),
                                    imageOutput(ns("hm_plot"), width = "100%", height = "800px")
                                )
                            )
                        )
                    ),
                    nav_panel("Results Table",
                        DTOutput(ns("results_table"))
                    ),
                    nav_panel("Common DEGs",
                        fluidRow(
                            column(3,
                                card(
                                    card_header("Parameters"),
                                    radioButtons(ns("cmp_type"), "Plot type",
                                        choices  = c("Venn Diagram"       = "venn",
                                                     "UpSet Plot"         = "upset",
                                                     "Log2FC Correlation" = "heatmap"),
                                        selected = "upset"),
                                    hr(),
                                    strong("Contrasts to compare"),
                                    uiOutput(ns("cmp_contrast_check_ui")),
                                    hr(),
                                    numericInput(ns("cmp_padj"), "adj. p-value \u2264",
                                        value = 0.05, min = 0, max = 1, step = 0.01),
                                    numericInput(ns("cmp_lfc"), "|log2FC| \u2265",
                                        value = 0.5, min = 0, max = NA, step = 0.1),
                                    hr(),
                                    actionButton(ns("cmp_apply"), "Apply Changes",
                                        class = "btn-primary w-100",
                                        icon  = icon("play")),
                                    downloadButton(ns("cmp_download_up"), "Download UP image",
                                        class = "btn-outline-secondary w-100 mt-2"),
                                    downloadButton(ns("cmp_download_down"), "Download DOWN image",
                                        class = "btn-outline-secondary w-100 mt-2"),
                                    uiOutput(ns("cmp_venn_warn"))
                                )
                            ),
                            column(9,
                                fluidRow(
                                    column(12,
                                        card(
                                            card_header(
                                                tags$strong("UP-regulated",
                                                    style = "color:#d73027")
                                            ),
                                            plotOutput(ns("cmp_up"),
                                                       height = "420px")
                                        )
                                    )
                                ),
                                fluidRow(
                                    column(12,
                                        card(
                                            card_header(
                                                tags$strong("DOWN-regulated",
                                                    style = "color:#4575b4")
                                            ),
                                            plotOutput(ns("cmp_down"),
                                                       height = "420px")
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}

mod_dge_server <- function(id, app_data) {
    moduleServer(id, function(input, output, session) {
        data <- reactive(app_data())

        observe({
            updateSelectInput(session, "contrast",
                              choices = names(data()$dge))
        })

        # Frozen params — only updates when Apply is clicked OR contrast changes.
        # This prevents every keystroke in numericInputs from triggering recomputation.
        params <- reactive({
            list(
                contrast   = input$contrast,
                padj_cut   = input$padj_cut,
                lfc_cut    = input$lfc_cut,
                n_label    = input$n_label,
                label_size = input$label_size
            )
        }) |> bindEvent(input$apply, input$contrast, ignoreNULL = TRUE, ignoreInit = FALSE)

        current_dge <- reactive({
            p <- params()
            req(p$contrast)
            df <- data()$dge[[p$contrast]]
            if (is.null(df)) return(data.frame())
            df %>%
                mutate(
                    direction = case_when(
                        !is.na(padj) & !is.na(log2FoldChange) &
                            padj < p$padj_cut & log2FoldChange >= p$lfc_cut  ~ "Up",
                        !is.na(padj) & !is.na(log2FoldChange) &
                            padj < p$padj_cut & log2FoldChange <= -p$lfc_cut ~ "Down",
                        TRUE ~ "NS"
                    )
                )
        })

        parse_genes <- function(x) {
            g <- unlist(strsplit(x %||% "", "[,;[:space:]]+"))
            unique(g[nzchar(g)])
        }

        collect_volcano_params <- function(contrast) {
            list(
                contrast         = contrast,
                padj_cut         = input$padj_cut,
                lfc_cut          = input$lfc_cut,
                n_label          = input$n_label,
                label_size       = input$label_size,
                must_genes       = input$must_genes,
                highlight_labels = isTRUE(input$highlight_labels),
                hide_top_n       = isTRUE(input$hide_top_n)
            )
        }

        build_volcano_plot <- function(df_raw, p) {
            req(nrow(df_raw) > 0)

            df <- df_raw %>%
                filter(!is.na(padj), !is.na(log2FoldChange)) %>%
                mutate(
                    padj_plot = pmax(padj, .Machine$double.xmin),
                    neglog10_padj = -log10(padj_plot),
                    direction = case_when(
                        padj < p$padj_cut & log2FoldChange >= p$lfc_cut  ~ "Up",
                        padj < p$padj_cut & log2FoldChange <= -p$lfc_cut ~ "Down",
                        TRUE ~ "NS"
                    ),
                    label_base = ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, gene_id)
                )

            if (nrow(df) == 0) {
                return(
                    ggplot() +
                        annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                                 label = "No DGE rows available for volcano plot") +
                        theme_void()
                )
            }

            n_lab  <- min(p$n_label, nrow(df))
            sig_df <- df %>% filter(direction != "NS")
            top_ids <- if (!isTRUE(p$hide_top_n) && nrow(sig_df) > 0 && n_lab > 0) {
                max_lfc <- max(abs(sig_df$log2FoldChange), na.rm = TRUE)
                max_nlp <- max(sig_df$neglog10_padj, na.rm = TRUE)
                sig_df %>%
                    mutate(score = sqrt(
                        (abs(log2FoldChange) / pmax(max_lfc, 1e-9))^2 +
                        (neglog10_padj / pmax(max_nlp, 1e-9))^2
                    )) %>%
                    arrange(desc(score)) %>%
                    head(n_lab) %>%
                    pull(gene_id)
            } else character(0)
            must <- parse_genes(p$must_genes)
            df$is_must <- toupper(df$label_base) %in% toupper(must) |
                toupper(df$gene_id) %in% toupper(must)
            highlight <- isTRUE(p$highlight_labels)

            # Highlighted genes get their own yellow layer, so skip them in the regular one
            regular_ids <- if (highlight) top_ids else union(top_ids, df$gene_id[df$is_must])
            df$label <- ifelse(df$gene_id %in% regular_ids & !(highlight & df$is_must),
                               df$label_base, NA_character_)
            hl_df <- df[df$is_must, , drop = FALSE]

            dir_col <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "#bbbbbb")

            gg <- ggplot(df, aes(x = log2FoldChange, y = neglog10_padj, colour = direction)) +
                geom_point(size = 0.8, alpha = 0.7) +
                ggrepel::geom_text_repel(aes(label = label),
                                         size = p$label_size / 3,
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.2,
                                         seed = 42,
                                         na.rm = TRUE)

            if (highlight && nrow(hl_df) > 0) {
                gg <- gg +
                    geom_point(data = hl_df,
                               aes(x = log2FoldChange, y = neglog10_padj),
                               inherit.aes = FALSE,
                               shape = 21, size = 2.8, stroke = 0.8,
                               fill = "#FFD700", colour = "grey20") +
                    ggrepel::geom_text_repel(data = hl_df,
                                             aes(x = log2FoldChange, y = neglog10_padj,
                                                 label = label_base),
                                             inherit.aes = FALSE,
                                             colour = "#FFD700",
                                             bg.color = "grey20", bg.r = 0.15,
                                             fontface = "bold",
                                             size = p$label_size / 3,
                                             max.overlaps = Inf,
                                             box.padding = 0.35,
                                             point.padding = 0.2,
                                             seed = 42)
            }

            gg +
                geom_vline(xintercept = c(-p$lfc_cut, p$lfc_cut), linetype = "dashed", colour = "grey40") +
                geom_hline(yintercept = -log10(pmax(p$padj_cut, .Machine$double.xmin)), linetype = "dashed", colour = "grey40") +
                scale_colour_manual(values = dir_col) +
                labs(title = paste0("Volcano: ", p$contrast),
                     x = "Shrunken log2 fold change",
                     y = "-log10(padj)",
                     colour = NULL) +
                theme_bw(base_size = 12) +
                theme(plot.margin = ggplot2::margin(t = 20, r = 10, b = 5, l = 5, unit = "pt")) +
                coord_cartesian(clip = "off")
        }

        # Keep the fast workflow: changing contrast shows pre-rendered pipeline PNG.
        # Recompute only when Apply is clicked for the current contrast.
        applied_volcano_params <- reactiveVal(NULL)

        observeEvent(input$contrast, {
            applied_volcano_params(NULL)
        }, ignoreInit = TRUE)

        observeEvent(input$apply, {
            applied_volcano_params(collect_volcano_params(input$contrast))

            must <- parse_genes(input$must_genes)
            df_raw <- data()$dge[[input$contrast]]
            if (length(must) > 0 && !is.null(df_raw)) {
                known <- toupper(c(df_raw$gene_name, df_raw$gene_id))
                missing <- must[!toupper(must) %in% known]
                if (length(missing) > 0) {
                    showNotification(paste("Genes not found:", paste(missing, collapse = ", ")),
                                     type = "warning")
                }
            }
        })

        output$volcano_plot <- renderImage({
            contrast <- input$contrast
            req(contrast)

            p_applied <- applied_volcano_params()
            apply_active <- !is.null(p_applied) && identical(p_applied$contrast, contrast)

            # 1) Default view: precomputed volcano from pipeline (instant, no recompute)
            if (!apply_active) {
                png_path <- get_dge_volcano_png(RESULTS_DIR, contrast)
                if (!is.null(png_path)) {
                    return(list(
                        src = png_path,
                        contentType = "image/png",
                        style = "width:100%; height:700px; object-fit:contain; background:white;"
                    ))
                }
            }

            # 2) Apply clicked (or missing precomputed PNG): build custom volcano
            df_raw <- data()$dge[[contrast]]
            if (is.null(df_raw) || nrow(df_raw) == 0) {
                tmpfile <- tempfile(fileext = ".png")
                p_empty <- ggplot() +
                    annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                             label = "No DGE results available for this contrast") +
                    theme_void()
                ggplot2::ggsave(tmpfile, p_empty, width = 8, height = 8, dpi = 150, bg = "white")
                return(list(
                    src = tmpfile,
                    contentType = "image/png",
                    style = "width:100%; height:700px; object-fit:contain; background:white;"
                ))
            }

            p_use <- if (apply_active) p_applied else collect_volcano_params(contrast)

            tmpfile <- tempfile(fileext = ".png")
            p_vol <- build_volcano_plot(df_raw, p_use)
            ggplot2::ggsave(tmpfile, p_vol, width = 8, height = 8, dpi = 150, bg = "white")
            list(
                src = tmpfile,
                contentType = "image/png",
                style = "width:100%; height:700px; object-fit:contain; background:white;"
            )
        }, deleteFile = FALSE)

        # Build the currently applied volcano plot; shared by the PNG/PDF exports.
        current_volcano_plot <- function(contrast) {
            df_raw <- data()$dge[[contrast]]
            if (is.null(df_raw) || nrow(df_raw) == 0) return(NULL)

            p_applied <- applied_volcano_params()
            p_use <- if (!is.null(p_applied) && identical(p_applied$contrast, contrast)) {
                p_applied
            } else {
                collect_volcano_params(contrast)
            }

            build_volcano_plot(df_raw, p_use)
        }

        output$download_volcano <- downloadHandler(
            filename = function() {
                contrast <- input$contrast %||% "contrast"
                paste0(contrast, "_volcano_dashboard.png")
            },
            content = function(file) {
                contrast <- input$contrast
                req(contrast)
                p_vol <- current_volcano_plot(contrast)
                req(!is.null(p_vol))
                ggplot2::ggsave(file, p_vol, device = "png",
                                width = 8, height = 8, dpi = 300, bg = "white")
            }
        )

        output$download_volcano_pdf <- downloadHandler(
            filename = function() {
                contrast <- input$contrast %||% "contrast"
                paste0(contrast, "_volcano_dashboard.pdf")
            },
            content = function(file) {
                contrast <- input$contrast
                req(contrast)
                p_vol <- current_volcano_plot(contrast)
                req(!is.null(p_vol))
                ggplot2::ggsave(file, p_vol, device = "pdf",
                                width = 8, height = 8, bg = "white")
            }
        )

        # =====================================================================
        # HEATMAP — server logic (ComplexHeatmap + design-aware VST per group)
        # =====================================================================
        hm_pkgs_ok <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
                      requireNamespace("circlize", quietly = TRUE)

        collect_hm_params <- function(contrast) {
            list(
                contrast     = contrast,
                top_n        = input$hm_top_n %||% 50,
                lfc          = input$hm_lfc %||% 1,
                padj         = input$hm_padj %||% 0.05,
                rank         = input$hm_rank %||% "logfc",
                mode         = input$hm_mode %||% "abs",
                names        = isTRUE(input$hm_names),
                cluster_rows = isTRUE(input$hm_cluster_rows),
                cluster_cols = isTRUE(input$hm_cluster_cols),
                dend         = isTRUE(input$hm_dend),
                supervise    = isTRUE(input$hm_supervise)
            )
        }

        # Frozen params: default (pre-rendered) view until Apply is pressed.
        applied_hm_params <- reactiveVal(NULL)
        observeEvent(input$contrast, { applied_hm_params(NULL) }, ignoreInit = TRUE)
        observeEvent(input$hm_apply, {
            applied_hm_params(collect_hm_params(input$contrast))
        })

        # Select top genes: significant (padj + |log2FC| cutoffs), ranked by the
        # chosen metric; "abs" = top N overall, "balanced" = N/2 up + N/2 down.
        select_heatmap_genes <- function(df, p) {
            df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), , drop = FALSE]
            df <- df[df$padj < p$padj & abs(df$log2FoldChange) >= p$lfc, , drop = FALSE]
            if (nrow(df) == 0) {
                return(list(genes = character(0),
                            note = "No significant genes for the current cutoffs."))
            }
            metric_of <- function(sub) {
                if (identical(p$rank, "padj")) sub$padj else abs(sub$log2FoldChange)
            }
            if (identical(p$mode, "balanced")) {
                n_half <- ceiling(p$top_n / 2)
                up <- df[df$log2FoldChange > 0, , drop = FALSE]
                dn <- df[df$log2FoldChange < 0, , drop = FALSE]
                up <- up[order(metric_of(up)), , drop = FALSE]
                dn <- dn[order(metric_of(dn)), , drop = FALSE]
                genes <- c(head(up$gene_id, n_half), head(dn$gene_id, n_half))
            } else {
                inc <- identical(p$rank, "padj")  # padj: ascending; |log2FC|: descending
                ord <- order(metric_of(df), decreasing = !inc)
                genes <- head(df$gene_id[ord], p$top_n)
            }
            list(genes = unique(genes), note = NULL)
        }

        # Build the ComplexHeatmap object, or a message-only placeholder.
        build_heatmap_obj <- function(contrast, p) {
            df <- data()$dge[[contrast]]
            if (is.null(df) || nrow(df) == 0) {
                return(list(ht = NULL, msg = "No DGE results for this contrast."))
            }
            sel <- select_heatmap_genes(df, p)
            if (length(sel$genes) < 2) {
                return(list(ht = NULL, msg = sel$note %||% "Too few genes to plot."))
            }

            # Design-aware VST for this contrast's norm_group (matches the pipeline).
            grp      <- get_contrast_norm_group(data()$contrasts_meta, contrast)
            vst_df   <- data()$norm_vst[[grp]]
            src_note <- NULL
            if (is.null(vst_df)) {
                vst_df   <- data()$vst
                src_note <- paste0("Design-aware VST for group '", grp,
                                   "' not found; using global blind VST.")
            }
            if (is.null(vst_df) || nrow(vst_df) == 0) {
                return(list(ht = NULL, msg = "No VST matrix available."))
            }

            gene_name_map <- setNames(vst_df$gene_name, vst_df$gene_id)
            vst_mat <- as.matrix(vst_df[, !colnames(vst_df) %in% c("gene_id", "gene_name"), drop = FALSE])
            rownames(vst_mat) <- vst_df$gene_id

            ss <- data()$samplesheet
            cm <- data()$contrasts_meta
            if (!is.null(cm) && nrow(cm) > 0 && contrast %in% cm$contrast_id) {
                row          <- cm[match(contrast, cm$contrast_id), ]
                contrast_var <- row$variable
                ref <- row$reference
                trt <- row$treatment
            } else {
                contrast_var <- "condition"
                ref <- trt <- NULL
            }
            if (is.null(ss) || nrow(ss) == 0) {
                return(list(ht = NULL, msg = "No samplesheet available."))
            }
            cond_col <- if (!is.null(contrast_var) && contrast_var %in% colnames(ss)) {
                ss[[contrast_var]]
            } else if ("condition" %in% colnames(ss)) {
                ss$condition
            } else {
                NULL
            }
            if (is.null(cond_col)) {
                return(list(ht = NULL, msg = "Samplesheet has no condition column."))
            }
            keep_ss <- if (!is.null(ref) && !is.null(trt)) cond_col %in% c(ref, trt) else rep(TRUE, nrow(ss))
            samples <- intersect(ss$sample[keep_ss], colnames(vst_mat))
            if (length(samples) < 2) {
                return(list(ht = NULL, msg = "Fewer than 2 samples for this contrast."))
            }

            mat  <- vst_mat[intersect(sel$genes, rownames(vst_mat)), samples, drop = FALSE]
            keep <- apply(mat, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
            mat  <- mat[keep, , drop = FALSE]
            if (nrow(mat) < 2) {
                return(list(ht = NULL, msg = "Too few variable genes for the selected cutoffs."))
            }
            mat <- mat - rowMeans(mat)  # center per row (matches pipeline)

            labs <- gene_name_map[rownames(mat)]
            miss <- is.na(labs) | !nzchar(labs)
            labs[miss] <- rownames(mat)[miss]
            rownames(mat) <- make.unique(labs)

            sample_cond <- cond_col[match(samples, ss$sample)]
            cond_vec    <- if (!is.null(ref) && !is.null(trt)) {
                factor(sample_cond, levels = unique(c(ref, trt)))
            } else {
                factor(sample_cond)
            }
            cond_lvls <- levels(cond_vec)
            n_pal     <- max(3, min(length(cond_lvls), 9))
            cond_cols <- setNames(
                colorRampPalette(brewer.pal(n_pal, "Set1"))(length(cond_lvls)),
                cond_lvls)

            lim     <- max(stats::quantile(abs(mat), 0.98, na.rm = TRUE), 1e-6)
            col_fun <- circlize::colorRamp2(c(-lim, 0, lim),
                                            c("#4575b4", "white", "#d73027"))

            cluster_cols  <- isTRUE(p$cluster_cols) && !isTRUE(p$supervise)
            show_col_dend <- cluster_cols && isTRUE(p$dend)
            show_row_dend <- isTRUE(p$cluster_rows) && isTRUE(p$dend)

            ht <- ComplexHeatmap::Heatmap(
                mat,
                name              = "Centered VST",
                col               = col_fun,
                top_annotation    = ComplexHeatmap::HeatmapAnnotation(
                    condition          = cond_vec,
                    col                = list(condition = cond_cols),
                    annotation_name_gp = grid::gpar(fontsize = 8)),
                show_row_names    = isTRUE(p$names),
                show_column_names = TRUE,
                row_names_gp      = grid::gpar(fontsize = 6),
                column_names_gp   = grid::gpar(fontsize = 8),
                cluster_rows      = isTRUE(p$cluster_rows),
                cluster_columns   = cluster_cols,
                show_row_dend     = show_row_dend,
                show_column_dend  = show_col_dend,
                column_split      = if (isTRUE(p$supervise)) cond_vec else NULL,
                column_title      = paste0("Top ", nrow(mat), " genes: ", contrast),
                column_title_gp   = grid::gpar(fontsize = 12, fontface = "bold"),
                heatmap_legend_param = list(title_gp  = grid::gpar(fontsize = 9),
                                            labels_gp = grid::gpar(fontsize = 8))
            )
            list(ht = ht, msg = src_note)
        }

        draw_heatmap_obj <- function(obj) {
            if (is.null(obj$ht)) {
                grid::grid.newpage()
                grid::grid.text(obj$msg %||% "Heatmap unavailable.",
                                gp = grid::gpar(col = "grey45", fontsize = 12))
            } else {
                ComplexHeatmap::draw(obj$ht)
            }
        }

        # Render a heatmap to a raster file at a given resolution.
        render_heatmap_file <- function(file, contrast, p, res = 150) {
            png(file, width = 8 * res, height = 8 * res, res = res)
            on.exit(dev.off(), add = TRUE)
            if (!hm_pkgs_ok) {
                grid::grid.newpage()
                grid::grid.text("ComplexHeatmap unavailable.\nPlease rebuild the container.",
                                gp = grid::gpar(col = "firebrick", fontsize = 12))
            } else {
                draw_heatmap_obj(build_heatmap_obj(contrast, p))
            }
        }

        current_hm_params <- function(contrast) {
            p_applied <- applied_hm_params()
            if (!is.null(p_applied) && identical(p_applied$contrast, contrast)) {
                p_applied
            } else {
                collect_hm_params(contrast)
            }
        }

        output$hm_plot <- renderImage({
            contrast <- input$contrast
            req(contrast)

            p_applied    <- applied_hm_params()
            apply_active <- !is.null(p_applied) && identical(p_applied$contrast, contrast)

            # 1) Default view: pipeline pre-rendered top-50 heatmap (instant, no compute)
            if (!apply_active) {
                png_path <- get_dge_heatmap_png(RESULTS_DIR, contrast, 50)
                if (!is.null(png_path)) {
                    return(list(
                        src         = png_path,
                        contentType = "image/png",
                        style       = "width:100%; height:800px; object-fit:contain; background:white;"))
                }
            }

            # 2) Apply clicked (or no pre-rendered PNG): recompute with ComplexHeatmap
            tmpfile <- tempfile(fileext = ".png")
            render_heatmap_file(tmpfile, contrast, current_hm_params(contrast), res = 120)
            list(
                src         = tmpfile,
                contentType = "image/png",
                style       = "width:100%; height:800px; object-fit:contain; background:white;")
        }, deleteFile = FALSE)

        output$hm_download <- downloadHandler(
            filename = function() {
                paste0(input$contrast %||% "contrast", "_heatmap_dashboard.png")
            },
            content = function(file) {
                contrast <- input$contrast
                req(contrast)
                render_heatmap_file(file, contrast, current_hm_params(contrast), res = 300)
            }
        )

        output$hm_download_pdf <- downloadHandler(
            filename = function() {
                paste0(input$contrast %||% "contrast", "_heatmap_dashboard.pdf")
            },
            content = function(file) {
                contrast <- input$contrast
                req(contrast)
                pdf(file, width = 8, height = 8)
                on.exit(dev.off(), add = TRUE)
                if (!hm_pkgs_ok) {
                    grid::grid.newpage()
                    grid::grid.text("ComplexHeatmap unavailable.\nPlease rebuild the container.",
                                    gp = grid::gpar(col = "firebrick", fontsize = 12))
                } else {
                    draw_heatmap_obj(build_heatmap_obj(contrast, current_hm_params(contrast)))
                }
            }
        )

        # =====================================================================
        # COMMON DEGs — server logic
        # Requires container >= 1.0.4 (ggVennDiagram, ComplexUpset,
        # ComplexHeatmap, circlize).
        # =====================================================================
        cmp_pkgs_ok <- all(sapply(
            c("ggVennDiagram", "ComplexUpset", "ComplexHeatmap", "circlize"),
            requireNamespace, quietly = TRUE
        ))

        # ---- Populate contrast checkboxes -----------------------------------
        output$cmp_contrast_check_ui <- renderUI({
            ns  <- session$ns
            cids <- names(data()$dge)
            if (length(cids) == 0) {
                return(p("No DGE results available.", style = "color:grey;font-size:0.9em;"))
            }
            checkboxGroupInput(ns("cmp_contrasts"), label = NULL,
                               choices = cids, selected = cids)
        })

        # ---- Venn safety warning: shown when >5 contrasts + venn selected ---
        output$cmp_venn_warn <- renderUI({
            req(input$cmp_type == "venn")
            n <- length(input$cmp_contrasts)
            if (!is.null(n) && n > 5) {
                div(class = "alert alert-warning mt-2",
                    style = "font-size:0.83em; padding:8px;",
                    tags$strong("\u26a0 Too many contrasts for a Venn Diagram (max 5)."),
                    tags$br(),
                    "Please select \u22645 contrasts or switch to UpSet Plot.")
            }
        })

        # ---- Main computation gated behind "Apply Changes" button -----------
        # ignoreNULL = FALSE: compute on app startup with default selections.
        cmp_sets <- eventReactive(input$cmp_apply, {
            req(data()$dge)
            selected <- input$cmp_contrasts
            validate(
                need(length(selected) >= 2,
                     "Please select at least 2 contrasts to compare.")
            )
            dge_list  <- data()$dge
            padj_thr  <- input$cmp_padj %||% 0.05
            lfc_thr   <- input$cmp_lfc  %||% 0.5

            get_sets <- function(direction) {
                lapply(setNames(selected, selected), function(cid) {
                    df <- dge_list[[cid]]
                    if (is.null(df) || nrow(df) == 0) return(character(0))
                    keep_padj <- !is.na(df$padj) & df$padj < padj_thr
                    keep_lfc  <- if (direction == "up") df$log2FoldChange >  lfc_thr
                                 else                   df$log2FoldChange < -lfc_thr
                    df$gene_id[keep_padj & keep_lfc]
                })
            }

            list(
                up       = get_sets("up"),
                down     = get_sets("down"),
                dge_full = dge_list[selected]
            )
        }, ignoreNULL = FALSE)

        # ---- Drawing helper (pure side-effect; works with renderPlot) -------
        draw_common_plot <- function(sets, dge_full, direction, plot_type) {
            colour <- if (direction == "up") "#d73027" else "#4575b4"
            label  <- if (direction == "up") "UP-regulated" else "DOWN-regulated"

            # Package guard
            if (!cmp_pkgs_ok) {
                grid::grid.text(
                    "Package(s) unavailable.\nPlease rebuild the container to version \u2265 1.0.5.",
                    gp = grid::gpar(fontsize = 12, col = "firebrick"))
                return(invisible(NULL))
            }

            # Minimum contrasts guard
            if (length(sets) < 2) {
                grid::grid.text("Select at least 2 contrasts.",
                                gp = grid::gpar(fontsize = 13, col = "grey50"))
                return(invisible(NULL))
            }

            if (plot_type == "venn") {
                # Enforce 5-set limit
                if (length(sets) > 5) {
                    grid::grid.text(
                        paste0("Too many contrasts for a Venn Diagram (max 5).\n",
                               "Please select \u22645 or switch to UpSet Plot."),
                        gp = grid::gpar(fontsize = 12, col = "firebrick"))
                    return(invisible(NULL))
                }
                if (all(sapply(sets, length) == 0)) {
                    grid::grid.text(paste0("No significant ", label, " DEGs found."),
                                    gp = grid::gpar(fontsize = 12, col = "grey50"))
                    return(invisible(NULL))
                }
                p <- ggVennDiagram::ggVennDiagram(
                        sets,
                        label     = "both",
                        set_color = rep("black", length(sets))) +
                    ggplot2::scale_fill_gradient(low = "white", high = colour) +
                    ggplot2::scale_color_manual(
                        values = rep("black", length(sets))) +
                    ggplot2::ggtitle(paste0(label, " DEGs")) +
                    ggplot2::theme(
                        legend.position = "none",
                        text = ggplot2::element_text(color = "black"),
                        plot.title = ggplot2::element_text(
                            hjust = 0.5, face = "bold", size = 13),
                        plot.margin = ggplot2::margin(16, 36, 16, 36, unit = "pt")) +
                    ggplot2::coord_cartesian(clip = "off")
                print(p)

            } else if (plot_type == "upset") {
                all_genes <- unique(unlist(sets))
                if (length(all_genes) == 0) {
                    grid::grid.text(paste0("No significant ", label, " DEGs found."),
                                    gp = grid::gpar(fontsize = 12, col = "grey50"))
                    return(invisible(NULL))
                }
                mat    <- as.data.frame(
                    sapply(names(sets), function(cid) all_genes %in% sets[[cid]]),
                    stringsAsFactors = FALSE)
                mat[]  <- lapply(mat, as.logical)
                p <- ComplexUpset::upset(
                    mat,
                    intersect    = names(sets),
                    name         = "Contrast",
                    set_sizes    = ComplexUpset::upset_set_size(
                        geom = ggplot2::geom_bar(fill = colour)),
                    width_ratio  = 0.2,
                    sort_intersections_by = "degree") +
                    ggplot2::ggtitle(paste0(label, " DEGs")) +
                    ggplot2::theme(
                        text = ggplot2::element_text(color = "black"),
                        axis.text = ggplot2::element_text(color = "black"),
                        axis.title = ggplot2::element_text(color = "black"),
                        plot.title = ggplot2::element_text(
                            hjust = 0.5, face = "bold", size = 13))
                print(p)

            } else if (plot_type == "heatmap") {
                sig_genes <- unique(unlist(sets))
                if (length(sig_genes) < 2) {
                    grid::grid.text(
                        paste0("Insufficient ", label, " DEGs for correlation heatmap."),
                        gp = grid::gpar(fontsize = 12, col = "grey50"))
                    return(invisible(NULL))
                }

                # Build log2FC matrix: rows = sig_genes, cols = contrasts
                lfc_mat <- sapply(names(dge_full), function(cid) {
                    df  <- dge_full[[cid]]
                    vec <- setNames(df$log2FoldChange, df$gene_id)
                    vec[sig_genes]
                })
                rownames(lfc_mat) <- sig_genes

                # Retain genes observed in >= 2 contrasts
                n_present <- rowSums(!is.na(lfc_mat))
                lfc_mat   <- lfc_mat[n_present >= 2, , drop = FALSE]

                if (nrow(lfc_mat) < 2 || ncol(lfc_mat) < 2) {
                    grid::grid.text(
                        "Insufficient gene overlap for heatmap.",
                        gp = grid::gpar(fontsize = 12, col = "grey50"))
                    return(invisible(NULL))
                }

                cor_mat <- cor(lfc_mat, method = "spearman",
                               use = "pairwise.complete.obs")
                cor_mat[is.nan(cor_mat)] <- NA

                col_fun <- circlize::colorRamp2(
                    c(-1, 0, 1), c("#4575b4", "white", "#d73027"))

                ht <- ComplexHeatmap::Heatmap(
                    cor_mat,
                    name              = "r",
                    col               = col_fun,
                    column_title      = paste0(
                        label, " DEGs - Spearman correlation",
                        "\n(n = ", nrow(lfc_mat), " genes)"),
                    column_title_gp   = grid::gpar(fontsize = 11, fontface = "bold"),
                    row_names_gp      = grid::gpar(fontsize = 9),
                    column_names_gp   = grid::gpar(fontsize = 9),
                    rect_gp           = grid::gpar(col = "white", lwd = 1),
                    na_col            = "#cccccc",
                    cluster_rows      = TRUE,
                    cluster_columns   = TRUE,
                    cell_fun          = function(j, i, x, y, width, height, fill) {
                        if (!is.na(cor_mat[i, j]))
                            grid::grid.text(
                                sprintf("%.2f", cor_mat[i, j]), x, y,
                                gp = grid::gpar(fontsize = 8, col = "black"))
                    }
                )
                ComplexHeatmap::draw(ht)
            }
        }

        # ---- Render UP plot -------------------------------------------------
        output$cmp_up <- renderPlot({
            d <- cmp_sets()
            req(d)
            draw_common_plot(d$up, d$dge_full, "up", input$cmp_type %||% "upset")
        })

        # ---- Render DOWN plot -----------------------------------------------
        output$cmp_down <- renderPlot({
            d <- cmp_sets()
            req(d)
            draw_common_plot(d$down, d$dge_full, "down", input$cmp_type %||% "upset")
        })

        output$cmp_download_up <- downloadHandler(
            filename = function() {
                paste0("common_degs_up_", input$cmp_type %||% "plot", ".png")
            },
            content = function(file) {
                d <- cmp_sets()
                req(d)
                png(file, width = 1400, height = 900, res = 180)
                on.exit(dev.off(), add = TRUE)
                draw_common_plot(d$up, d$dge_full, "up", input$cmp_type %||% "upset")
            }
        )

        output$cmp_download_down <- downloadHandler(
            filename = function() {
                paste0("common_degs_down_", input$cmp_type %||% "plot", ".png")
            },
            content = function(file) {
                d <- cmp_sets()
                req(d)
                png(file, width = 1400, height = 900, res = 180)
                on.exit(dev.off(), add = TRUE)
                draw_common_plot(d$down, d$dge_full, "down", input$cmp_type %||% "upset")
            }
        )

        output$results_table <- renderDT({
            df <- current_dge()
            if (nrow(df) == 0) return(datatable(data.frame()))
            df_show <- df %>%
                mutate(
                    log2FoldChange = round(log2FoldChange, 3),
                    baseMean       = round(baseMean, 1),
                    pvalue         = signif(pvalue, 3),
                    padj           = signif(padj, 3)
                )
            df_show <- df_show[, c("gene_id", "gene_name", "baseMean", "log2FoldChange",
                                   "pvalue", "padj", "direction"), drop = FALSE]

            datatable(df_show,
                      filter   = "top",
                      options  = list(pageLength = 20, scrollX = TRUE),
                      rownames = FALSE,
                      class    = "compact stripe") %>%
                formatStyle("direction",
                            backgroundColor = styleEqual(
                                c("Up", "Down", "NS"),
                                c("#fdd0c8", "#c8d8f8", "white")
                            ))
        })
    })
}
