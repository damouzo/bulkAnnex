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
                    uiOutput(ns("gene_selector_ui")),
                    selectizeInput(ns("selected_contrasts"), "Contrasts shown",
                                   choices = NULL, multiple = TRUE,
                                   options = list(plugins = list("remove_button"),
                                                  placeholder = "All contrasts"))
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Expression across conditions",
                        plotlyOutput(ns("expression_plot"), height = "450px"),
                        div(class = "mt-2",
                            downloadButton(ns("dl_expr_png"), "Download PNG",
                                           class = "btn-outline-secondary btn-sm"),
                            downloadButton(ns("dl_expr_pdf"), "Download PDF",
                                           class = "btn-outline-secondary btn-sm")
                        )
                    ),
                    nav_panel("DGE across contrasts",
                        plotlyOutput(ns("dge_lfc_plot"), height = "400px"),
                        div(class = "mt-2",
                            downloadButton(ns("dl_dge_png"), "Download PNG",
                                           class = "btn-outline-secondary btn-sm"),
                            downloadButton(ns("dl_dge_pdf"), "Download PDF",
                                           class = "btn-outline-secondary btn-sm")
                        ),
                        DTOutput(ns("dge_results_table"))
                    ),
                    nav_panel("GSEA membership",
                        DTOutput(ns("gsea_membership_table"))
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

        # ---- Contrast selector --------------------------------------------------
        # Populated from the available DGE contrasts. Defaults to all selected so
        # the sub-tabs behave exactly as before until the user narrows the set.
        observe({
            all_c <- names(data()$dge)
            updateSelectizeInput(session, "selected_contrasts",
                                 choices = all_c, selected = all_c)
        })

        # Active contrasts. An empty selection yields character(0) so the
        # sub-tabs can show a "no contrasts selected" notice.
        active_contrasts <- reactive({
            all_c <- names(data()$dge)
            sel   <- input$selected_contrasts
            if (is.null(sel)) all_c else intersect(all_c, sel)
        })

        # Conditions implied by the selected contrasts (via contrasts_meta).
        # NULL means "cannot filter" (legacy runs without contrasts_meta) and the
        # expression plot shows every condition.
        active_conditions <- reactive({
            cm  <- data()$contrasts_meta
            sel <- active_contrasts()
            if (is.null(cm) || nrow(cm) == 0) return(NULL)
            rows <- cm[cm$contrast_id %in% sel, , drop = FALSE]
            if (nrow(rows) == 0) return(NULL)
            unique(c(rows$reference, rows$treatment))
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

        # Expression data filtered by the active contrasts' conditions.
        gene_expr_filtered <- reactive({
            df <- gene_expr()
            if (is.null(df)) return(NULL)
            conds <- active_conditions()
            if (!is.null(conds)) df <- df[df$condition %in% conds, , drop = FALSE]
            df
        })

        selected_gene_name <- reactive({
            req(input$selected_gene)
            data()$vst$gene_name[data()$vst$gene_id == input$selected_gene][1]
        })

        # Safe label for download filenames — never errors when no gene is selected.
        gene_label <- function() {
            g <- tryCatch(selected_gene_name(), error = function(e) NULL)
            if (is.null(g) || length(g) == 0 || is.na(g) || !nzchar(g)) "gene" else g
        }

        output$expression_plot <- renderPlotly({
            if (length(active_contrasts()) == 0)
                return(plot_ly() %>% layout(title = "No contrasts selected"))
            df <- gene_expr_filtered()
            if (is.null(df)) return(plot_ly() %>% layout(title = "No gene selected"))
            if (nrow(df) == 0)
                return(plot_ly() %>% layout(title = "No conditions for the selected contrasts"))

            gene_name <- selected_gene_name()
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

        # ggplot equivalent of the expression boxplot for PNG/PDF export.
        build_expression_ggplot <- function(df, gene_name) {
            ggplot2::ggplot(df, ggplot2::aes(x = condition, y = vst_expr)) +
                ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
                ggplot2::geom_jitter(width = 0.15, height = 0, size = 2,
                                     alpha = 0.8, colour = "#6F42C1") +
                ggplot2::labs(title = paste0(gene_name, " \u2014 VST expression"),
                              x = "Condition", y = "VST expression") +
                ggplot2::theme_bw(base_size = 11)
        }

        output$dl_expr_png <- downloadHandler(
            filename = function() paste0(gene_label(), "_expression.png"),
            content  = function(file) {
                df <- gene_expr_filtered()
                req(df, nrow(df) > 0)
                ggplot2::ggsave(file, build_expression_ggplot(df, selected_gene_name()),
                                device = "png", width = 8, height = 6, dpi = 300, bg = "white")
            }
        )

        output$dl_expr_pdf <- downloadHandler(
            filename = function() paste0(gene_label(), "_expression.pdf"),
            content  = function(file) {
                df <- gene_expr_filtered()
                req(df, nrow(df) > 0)
                ggplot2::ggsave(file, build_expression_ggplot(df, selected_gene_name()),
                                device = "pdf", width = 8, height = 6, bg = "white")
            }
        )

        # Per-contrast log2FC for the selected gene (shared by plot and exports).
        dge_lfc_data <- reactive({
            req(input$selected_gene)
            dge_list <- data()$dge
            cids     <- active_contrasts()
            if (length(dge_list) == 0 || length(cids) == 0) return(NULL)

            lfc_df <- do.call(rbind, lapply(cids, function(cid) {
                df  <- dge_list[[cid]]
                row <- df[df$gene_id == input$selected_gene, , drop = FALSE]
                if (nrow(row) == 0) return(NULL)
                data.frame(
                    contrast = cid,
                    lfc      = row$log2FoldChange[1],
                    padj     = row$padj[1],
                    stringsAsFactors = FALSE
                )
            }))

            if (is.null(lfc_df) || nrow(lfc_df) == 0) return(NULL)

            lfc_df$sig   <- lfc_df$padj < 0.05 & !is.na(lfc_df$padj)
            lfc_df$color <- ifelse(lfc_df$sig, ifelse(lfc_df$lfc > 0, "#d73027", "#4575b4"), "#bbbbbb")
            lfc_df
        })

        output$dge_lfc_plot <- renderPlotly({
            req(input$selected_gene)
            if (length(active_contrasts()) == 0)
                return(plot_ly() %>% layout(title = "No contrasts selected"))
            lfc_df <- dge_lfc_data()
            if (is.null(lfc_df))
                return(plot_ly() %>% layout(title = "Gene not found in DGE results"))

            gene_name <- selected_gene_name()

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

        # ggplot equivalent of the log2FC bar chart for PNG/PDF export.
        build_dge_lfc_ggplot <- function(lfc_df, gene_name) {
            lfc_df$contrast <- factor(lfc_df$contrast,
                                      levels = lfc_df$contrast[order(lfc_df$lfc)])
            ggplot2::ggplot(lfc_df, ggplot2::aes(x = lfc, y = contrast, fill = color)) +
                ggplot2::geom_col() +
                ggplot2::scale_fill_identity() +
                ggplot2::geom_vline(xintercept = 0, colour = "black") +
                ggplot2::labs(title = paste0(gene_name, " \u2014 log2FC across contrasts"),
                              x = "log2 Fold Change", y = NULL) +
                ggplot2::theme_bw(base_size = 11)
        }

        output$dl_dge_png <- downloadHandler(
            filename = function() paste0(gene_label(), "_log2FC.png"),
            content  = function(file) {
                lfc_df <- dge_lfc_data()
                req(lfc_df)
                ggplot2::ggsave(file, build_dge_lfc_ggplot(lfc_df, selected_gene_name()),
                                device = "png", width = 8, height = 6, dpi = 300, bg = "white")
            }
        )

        output$dl_dge_pdf <- downloadHandler(
            filename = function() paste0(gene_label(), "_log2FC.pdf"),
            content  = function(file) {
                lfc_df <- dge_lfc_data()
                req(lfc_df)
                ggplot2::ggsave(file, build_dge_lfc_ggplot(lfc_df, selected_gene_name()),
                                device = "pdf", width = 8, height = 6, bg = "white")
            }
        )

        output$dge_results_table <- renderDT({
            req(input$selected_gene)
            dge_list <- data()$dge
            cids     <- active_contrasts()
            all_rows <- do.call(rbind, lapply(cids, function(cid) {
                df  <- dge_list[[cid]]
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
            active    <- active_contrasts()
            rows <- lapply(names(gsea_data), function(db) {
                db_data <- gsea_data[[db]]
                cids    <- intersect(names(db_data), active)
                lapply(cids, function(cid) {
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

    })
}
