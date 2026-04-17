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

            expr_row <- vst[vst$gene_id == input$selected_gene, , drop = FALSE]
            sample_cols <- ss$sample
            expr_vals <- as.numeric(expr_row[1, sample_cols])

            df <- data.frame(
                sample    = sample_cols,
                vst_expr  = expr_vals,
                condition = ss$condition,
                stringsAsFactors = FALSE
            )
            if ("batch" %in% colnames(ss)) df$batch <- ss$batch
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

                    mask <- grepl(gene_name, df[[col]], ignore.case = TRUE) |
                            grepl(gene_id,   df[[col]], ignore.case = TRUE)
                    if (!any(mask, na.rm = TRUE)) return(NULL)

                    cols_exist <- intersect(c("Description", "NES", "padj"), colnames(df))
                    if (length(cols_exist) == 0) return(NULL)
                    df_hit <- df[mask, cols_exist, drop = FALSE]
                    cbind(database = db, contrast = cid, df_hit)
                })
            })

            all_rows <- do.call(rbind, Filter(Negate(is.null),
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
