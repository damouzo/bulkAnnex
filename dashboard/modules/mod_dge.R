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
                    actionButton(ns("apply"), "Apply settings",
                                 icon = icon("play"),
                                 class = "btn-primary w-100 mt-2"),
                    downloadButton(ns("download_volcano"), "Download volcano PNG",
                                   class = "btn-outline-secondary w-100 mt-2")
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Volcano",
                        imageOutput(ns("volcano_plot"), width = "100%", height = "700px")
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
            top_ids <- if (nrow(sig_df) > 0 && n_lab > 0) {
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
            df$label <- ifelse(df$gene_id %in% top_ids, df$label_base, NA_character_)

            dir_col <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "#bbbbbb")

            ggplot(df, aes(x = log2FoldChange, y = neglog10_padj, colour = direction)) +
                geom_point(size = 0.8, alpha = 0.7) +
                ggrepel::geom_text_repel(aes(label = label),
                                         size = p$label_size / 3,
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.2,
                                         seed = 42,
                                         na.rm = TRUE) +
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
            applied_volcano_params(list(
                contrast = input$contrast,
                padj_cut = input$padj_cut,
                lfc_cut = input$lfc_cut,
                n_label = input$n_label,
                label_size = input$label_size
            ))
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

            p_use <- if (apply_active) p_applied else list(
                contrast = contrast,
                padj_cut = input$padj_cut,
                lfc_cut = input$lfc_cut,
                n_label = input$n_label,
                label_size = input$label_size
            )

            tmpfile <- tempfile(fileext = ".png")
            p_vol <- build_volcano_plot(df_raw, p_use)
            ggplot2::ggsave(tmpfile, p_vol, width = 8, height = 8, dpi = 150, bg = "white")
            list(
                src = tmpfile,
                contentType = "image/png",
                style = "width:100%; height:700px; object-fit:contain; background:white;"
            )
        }, deleteFile = FALSE)

        output$download_volcano <- downloadHandler(
            filename = function() {
                contrast <- input$contrast %||% "contrast"
                paste0(contrast, "_volcano_dashboard.png")
            },
            content = function(file) {
                contrast <- input$contrast
                req(contrast)
                df_raw <- data()$dge[[contrast]]
                req(!is.null(df_raw), nrow(df_raw) > 0)

                p_applied <- applied_volcano_params()
                p_use <- if (!is.null(p_applied) && identical(p_applied$contrast, contrast)) {
                    p_applied
                } else {
                    list(
                        contrast = contrast,
                        padj_cut = input$padj_cut,
                        lfc_cut = input$lfc_cut,
                        n_label = input$n_label,
                        label_size = input$label_size
                    )
                }

                p_vol <- build_volcano_plot(df_raw, p_use)
                ggplot2::ggsave(file, p_vol, width = 8, height = 8, dpi = 300, bg = "white")
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
