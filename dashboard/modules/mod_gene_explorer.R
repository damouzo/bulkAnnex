# modules/mod_gene_explorer.R — Gene Explorer tab

mod_gene_explorer_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12, h3("Gene Explorer"))
        ),
        fluidRow(
            column(3,
                card(
                    card_header("Search"),
                    textInput(ns("gene_query"), "Gene name or Ensembl ID",
                              placeholder = "e.g. TP53, ENSG00000141510"),
                    actionButton(ns("search_btn"), "Search", class = "btn-primary"),
                    hr(),
                    uiOutput(ns("gene_selector_ui"))
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Expression across conditions",
                        plotlyOutput(ns("expression_plot"), height = "450px")
                    ),
                    nav_panel("DGE across contrasts",
                        plotlyOutput(ns("dge_lfc_plot"), height = "400px"),
                        DTOutput(ns("dge_results_table"))
                    ),
                    nav_panel("GSEA membership",
                        DTOutput(ns("gsea_membership_table"))
                    ),
                    # ---- Common DEGs sub-tab --------------------------------
                    nav_panel("Common DEGs",
                        fluidRow(
                            # Controls column
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
                                    sliderInput(ns("cmp_padj"), "adj. p-value \u2264",
                                        min = 0.001, max = 0.2, value = 0.05, step = 0.005),
                                    sliderInput(ns("cmp_lfc"), "|log2FC| \u2265",
                                        min = 0, max = 3, value = 0.5, step = 0.1),
                                    hr(),
                                    actionButton(ns("cmp_apply"), "Apply Changes",
                                        class = "btn-primary w-100",
                                        icon  = icon("play")),
                                    uiOutput(ns("cmp_venn_warn"))
                                )
                            ),
                            # Plot area: UP left, DOWN right
                            column(9,
                                fluidRow(
                                    column(6,
                                        card(
                                            card_header(
                                                tags$strong("UP-regulated",
                                                    style = "color:#d73027")
                                            ),
                                            plotOutput(ns("cmp_up"),
                                                       height = "460px")
                                        )
                                    ),
                                    column(6,
                                        card(
                                            card_header(
                                                tags$strong("DOWN-regulated",
                                                    style = "color:#4575b4")
                                            ),
                                            plotOutput(ns("cmp_down"),
                                                       height = "460px")
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

mod_gene_explorer_server <- function(id, app_data) {
    moduleServer(id, function(input, output, session) {
        data <- reactive(app_data())

        # Search results
        search_results <- eventReactive(input$search_btn, {
            req(input$gene_query, nchar(trimws(input$gene_query)) > 0)
            vst <- data()$vst
            if (nrow(vst) == 0) return(data.frame())
            q <- trimws(input$gene_query)
            mask <- grepl(q, vst$gene_id, ignore.case = TRUE) |
                    grepl(q, vst$gene_name, ignore.case = TRUE)
            vst[mask, c("gene_id", "gene_name"), drop = FALSE]
        })

        output$gene_selector_ui <- renderUI({
            ns <- session$ns
            # Before first search, search_results() throws silently via req() — show prompt
            res <- tryCatch(search_results(), error = function(e) NULL)
            if (is.null(res)) {
                return(p("Enter a gene name and click Search.", style = "color: grey; font-size: 0.9em;"))
            }
            if (nrow(res) == 0) {
                return(p("No genes found.", style = "color: grey; font-size: 0.9em;"))
            }
            choices <- setNames(res$gene_id,
                                paste0(res$gene_name, " (", res$gene_id, ")"))
            selectInput(ns("selected_gene"), "Select gene", choices = choices)
        })

        # Expression data for selected gene
        gene_expr <- reactive({
            req(input$selected_gene)
            vst <- data()$vst
            ss  <- data()$samplesheet
            if (nrow(vst) == 0 || !input$selected_gene %in% vst$gene_id) return(NULL)

            expr_row    <- vst[vst$gene_id == input$selected_gene, , drop = FALSE]
            # Only include samples present in both samplesheet and VST (multi-group safe)
            sample_cols <- intersect(ss$sample, colnames(vst))
            if (length(sample_cols) == 0) return(NULL)
            expr_vals <- as.numeric(expr_row[1, sample_cols])

            df <- data.frame(
                sample    = sample_cols,
                vst_expr  = expr_vals,
                condition = ss$condition[match(sample_cols, ss$sample)],
                stringsAsFactors = FALSE
            )
            if ("batch" %in% colnames(ss)) df$batch <- ss$batch[match(sample_cols, ss$sample)]
            df
        })

        output$expression_plot <- renderPlotly({
            df <- gene_expr()
            if (is.null(df)) return(plot_ly() %>% layout(title = "No gene selected"))

            gene_name <- data()$vst$gene_name[data()$vst$gene_id == input$selected_gene][1]
            title_str <- paste0(gene_name, " — VST expression")

            # Jittered dot + boxplot
            plot_ly(df,
                    x    = ~condition,
                    y    = ~vst_expr,
                    text = ~paste0("Sample: ", sample,
                                   "\nVST: ", round(vst_expr, 3),
                                   "\nCondition: ", condition),
                    hoverinfo = "text",
                    type      = "box",
                    boxpoints = "all",
                    jitter    = 0.4,
                    pointpos  = 0,
                    marker    = list(size = 6, opacity = 0.8)) %>%
                layout(
                    xaxis = list(title = "Condition"),
                    yaxis = list(title = "VST expression"),
                    title = title_str
                )
        })

        output$dge_lfc_plot <- renderPlotly({
            req(input$selected_gene)
            dge_list <- data()$dge
            if (length(dge_list) == 0) return(plot_ly() %>% layout(title = "No DGE data"))

            gene_name <- data()$vst$gene_name[data()$vst$gene_id == input$selected_gene][1]

            lfc_df <- do.call(rbind, lapply(names(dge_list), function(cid) {
                df <- dge_list[[cid]]
                row <- df[df$gene_id == input$selected_gene, , drop = FALSE]
                if (nrow(row) == 0) return(NULL)
                data.frame(
                    contrast = cid,
                    lfc      = row$log2FoldChange[1],
                    padj     = row$padj[1],
                    stringsAsFactors = FALSE
                )
            }))

            if (is.null(lfc_df) || nrow(lfc_df) == 0) {
                return(plot_ly() %>% layout(title = "Gene not found in DGE results"))
            }

            lfc_df$sig   <- lfc_df$padj < 0.05 & !is.na(lfc_df$padj)
            lfc_df$color <- ifelse(lfc_df$sig, ifelse(lfc_df$lfc > 0, "#d73027", "#4575b4"), "#bbbbbb")

            plot_ly(lfc_df,
                    x    = ~lfc,
                    y    = ~reorder(contrast, lfc),
                    text = ~paste0("padj: ", signif(padj, 3)),
                    hoverinfo  = "text",
                    type       = "bar",
                    orientation = "h",
                    marker     = list(color = ~color)) %>%
                layout(
                    xaxis = list(title = "log2 Fold Change"),
                    yaxis = list(title = NULL),
                    title = paste0(gene_name, " — log2FC across contrasts"),
                    shapes = list(list(type = "line", x0 = 0, x1 = 0,
                                       y0 = 0, y1 = 1, yref = "paper",
                                       line = list(color = "black")))
                )
        })

        output$dge_results_table <- renderDT({
            req(input$selected_gene)
            dge_list <- data()$dge
            all_rows <- do.call(rbind, lapply(names(dge_list), function(cid) {
                df <- dge_list[[cid]]
                row <- df[df$gene_id == input$selected_gene, , drop = FALSE]
                if (nrow(row) == 0) return(NULL)
                cbind(contrast = cid, row[, c("gene_name", "baseMean",
                                               "log2FoldChange", "pvalue", "padj")])
            }))
            if (is.null(all_rows)) return(datatable(data.frame()))
            all_rows$log2FoldChange <- round(all_rows$log2FoldChange, 3)
            all_rows$baseMean       <- round(all_rows$baseMean, 1)
            all_rows$pvalue         <- signif(all_rows$pvalue, 3)
            all_rows$padj           <- signif(all_rows$padj,   3)
            datatable(all_rows, rownames = FALSE, class = "compact stripe",
                      options = list(pageLength = 10))
        })

        output$gsea_membership_table <- renderDT({
            req(input$selected_gene)
            gene_id   <- input$selected_gene
            gene_name <- data()$vst$gene_name[data()$vst$gene_id == gene_id][1]

            gsea_data <- data()$gsea
            rows <- lapply(names(gsea_data), function(db) {
                db_data <- gsea_data[[db]]
                lapply(names(db_data), function(cid) {
                    df <- db_data[[cid]]
                    if (is.null(df) || nrow(df) == 0) return(NULL)

                    # Normalise Description column
                    if (!"Description" %in% colnames(df) && "pathway" %in% colnames(df))
                        df$Description <- df$pathway

                    # Normalise padj column
                    if (!"padj" %in% colnames(df) && "p.adjust" %in% colnames(df))
                        df$padj <- df$p.adjust
                    if (!"padj" %in% colnames(df) && "pval" %in% colnames(df))
                        df$padj <- p.adjust(df$pval, method = "BH")

                    # Find leading-edge column
                    le_col <- intersect(c("core_enrichment", "leadingEdge", "leading_edge"), colnames(df))
                    if (length(le_col) == 0) return(NULL)
                    col <- le_col[1]

                    # Determine search terms: always try gene_id (Ensembl) and gene_name (symbol).
                    # For KEGG/Reactome with old results (core_enrichment = Entrez IDs),
                    # also try to resolve the Entrez ID via AnnotationDbi if available.
                    search_terms <- c(gene_id, gene_name)
                    if (all(grepl("^[0-9/]+$", na.omit(df[[col]])))) {
                        # core_enrichment looks like Entrez IDs — try live mapping
                        entrez_id <- tryCatch({
                            if (requireNamespace("AnnotationDbi", quietly = TRUE)) {
                                org_pkg <- if (exists("org.Hs.eg.db", envir = .GlobalEnv, inherits = TRUE) ||
                                               "org.Hs.eg.db" %in% rownames(installed.packages()))
                                               get("org.Hs.eg.db",
                                                   envir = asNamespace("org.Hs.eg.db"))
                                           else NULL
                                if (!is.null(org_pkg))
                                    AnnotationDbi::mapIds(org_pkg, keys = gene_id,
                                                          column = "ENTREZID",
                                                          keytype = "ENSEMBL",
                                                          multiVals = "first")
                                else NA_character_
                            } else NA_character_
                        }, error = function(e) NA_character_)
                        if (!is.na(entrez_id)) search_terms <- c(search_terms, entrez_id)
                    }

                    mask <- Reduce(`|`, lapply(search_terms, function(term) {
                        grepl(term, df[[col]], ignore.case = TRUE, fixed = TRUE)
                    }))
                    if (!any(mask, na.rm = TRUE)) return(NULL)

                    cols_exist <- intersect(c("Description", "ontology", "NES", "padj"), colnames(df))
                    if (length(cols_exist) == 0) return(NULL)
                    df_hit <- df[mask, cols_exist, drop = FALSE]
                    # Ensure every row has the same fixed column set (NA for missing cols).
                    # This prevents rbind() failing when GO rows have 'ontology' but
                    # KEGG/Reactome rows do not.
                    for (col_std in c("Description", "ontology", "NES", "padj")) {
                        if (!col_std %in% colnames(df_hit))
                            df_hit[[col_std]] <- NA
                    }
                    df_hit <- df_hit[, c("Description", "ontology", "NES", "padj"), drop = FALSE]
                    data.frame(database = db, contrast = cid, df_hit,
                               stringsAsFactors = FALSE, check.names = FALSE)
                })
            })

            all_rows <- dplyr::bind_rows(Filter(Negate(is.null),
                                                unlist(rows, recursive = FALSE)))
            if (is.null(all_rows) || nrow(all_rows) == 0) {
                return(datatable(
                    data.frame(info = paste0(gene_name, " not found in any leading edge.")),
                    rownames = FALSE, options = list(dom = "t")
                ))
            }
            if ("NES"  %in% colnames(all_rows)) all_rows$NES  <- round(all_rows$NES, 3)
            if ("padj" %in% colnames(all_rows)) all_rows$padj <- signif(all_rows$padj, 3)
            datatable(all_rows, rownames = FALSE, class = "compact stripe",
                      options = list(pageLength = 15, scrollX = TRUE))
        })

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
                    "Package(s) unavailable.\nPlease rebuild the container to version \u2265 1.0.4.",
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
                        set_color = rep(colour, length(sets))) +
                    ggplot2::scale_fill_gradient(low = "white", high = colour) +
                    ggplot2::scale_color_manual(
                        values = rep(colour, length(sets))) +
                    ggplot2::ggtitle(paste0(label, " DEGs")) +
                    ggplot2::theme(
                        legend.position = "none",
                        plot.title = ggplot2::element_text(
                            hjust = 0.5, face = "bold", size = 13))
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
                        label, " DEGs \u2014 Spearman correlation",
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

    })
}
