# modules/mod_gsea.R — GSEA tab

mod_gsea_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12, h3("Gene Set Enrichment Analysis"))
        ),
        fluidRow(
            column(3,
                card(
                    card_header("Settings"),
                    selectInput(ns("contrast"),  "Contrast",  choices = NULL),
                    selectInput(ns("database"),  "Database",
                                choices = c("GO"       = "go",
                                            "KEGG"     = "kegg",
                                            "Hallmarks"= "hallmarks",
                                            "Reactome" = "reactome")),
                    # GO ontology sub-selector (only visible for GO)
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'go'", ns("database")),
                        selectInput(ns("go_ont"), "GO Ontology",
                                    choices = c("Biological Process" = "BP",
                                                "Molecular Function" = "MF",
                                                "Cellular Component" = "CC"),
                                    selected = "BP")
                    ),
                    sliderInput(ns("padj_cut"), "adj.p cutoff",
                                min = 0.001, max = 0.5, value = 0.05, step = 0.001),
                    numericInput(ns("n_terms"), "Top N terms",
                                 value = 20, min = 5, max = 50),
                    # Pathway selector — used by the GSEA Running Sum tab
                    selectInput(ns("pathway"), "Pathway (running sum)", choices = NULL)
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Dot Plot",
                        plotOutput(ns("dotplot"), height = "600px")
                    ),
                    nav_panel("Ridge Plot",
                        plotOutput(ns("ridgeplot"), height = "600px")
                    ),
                    nav_panel("Running Sum",
                        plotOutput(ns("gseaplot"), height = "500px")
                    ),
                    nav_panel("Results Table",
                        DTOutput(ns("results_table"))
                    )
                )
            )
        )
    )
}

mod_gsea_server <- function(id, app_data) {
    moduleServer(id, function(input, output, session) {
        data <- reactive(app_data())

        # Populate contrast selector
        observe({
            gsea_contrasts <- unique(unlist(lapply(data()$gsea, names)))
            choices <- if (length(gsea_contrasts) > 0) gsea_contrasts else names(data()$dge)
            updateSelectInput(session, "contrast", choices = choices)
        })

        # ---- Resolve the native gseaResult object (or NULL) -----------------
        # Returns:  gseaResult S4 object  — use enrichplot functions
        #           data.frame/data.table — use plotly fallback (Hallmarks / no RDS)
        #           NULL                  — no data available
        current_gse_obj <- reactive({
            req(input$contrast, input$database)
            objs <- data()$gsea_objects[[input$contrast]]
            if (is.null(objs)) return(NULL)
            if (input$database == "go") {
                ont_key <- paste0("go_", tolower(input$go_ont %||% "BP"))
                objs[[ont_key]]
            } else {
                objs[[input$database]]
            }
        })

        # ---- CSV-based data frame (used for Results Table + plotly fallback) -
        current_gsea_df <- reactive({
            req(input$contrast, input$database)
            # For GO, the CSV is the combined all-ontologies file
            db_data <- data()$gsea[[input$database]]
            if (is.null(db_data)) return(data.frame())
            df <- db_data[[input$contrast]]
            if (is.null(df) || nrow(df) == 0) return(data.frame())
            # Normalise column names for consistency
            if (!"padj" %in% colnames(df) && "pval" %in% colnames(df))
                df$padj <- p.adjust(df$pval, method = "BH")
            if (!"Description" %in% colnames(df) && "pathway" %in% colnames(df))
                df$Description <- df$pathway
            df
        })

        # When the selected contrast/db changes, populate pathway selector
        observe({
            obj <- current_gse_obj()
            choices <- character(0)
            if (!is.null(obj) && isS4(obj)) {
                df <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
                if (nrow(df) > 0 && "Description" %in% colnames(df))
                    choices <- df$Description[order(df$p.adjust)[seq_len(min(50, nrow(df)))]]
            } else {
                df <- current_gsea_df()
                if (nrow(df) > 0 && "Description" %in% colnames(df))
                    choices <- df$Description[order(df$padj)[seq_len(min(50, nrow(df)))]]
            }
            updateSelectInput(session, "pathway", choices = choices,
                              selected = if (length(choices) > 0) choices[1] else NULL)
        })

        # ---- Dot Plot -------------------------------------------------------
        output$dotplot <- renderPlot({
            obj <- current_gse_obj()
            if (!is.null(obj) && isS4(obj)) {
                # enrichplot::dotplot — uses the native gseaResult object
                tryCatch(
                    enrichplot::dotplot(obj, showCategory = input$n_terms,
                                        title = paste(input$database, "GSEA —", input$contrast)) +
                        theme_bw(base_size = 11),
                    error = function(e) {
                        ggplot() + annotate("text", x = 0.5, y = 0.5,
                                            label = paste("enrichplot error:", e$message)) +
                            theme_void()
                    }
                )
            } else {
                # Plotly-free ggplot fallback from CSV data
                df <- current_gsea_df()
                if (nrow(df) == 0)
                    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No results") + theme_void())
                df_plot <- df %>%
                    filter(padj < input$padj_cut) %>%
                    arrange(padj) %>%
                    head(input$n_terms) %>%
                    mutate(
                        Description = stringr::str_wrap(Description, 50),
                        gs_size = if ("size" %in% colnames(.)) size
                                  else if ("setSize" %in% colnames(.)) setSize
                                  else 50
                    )
                if (nrow(df_plot) == 0)
                    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No significant terms") + theme_void())
                ggplot(df_plot, aes(x = NES, y = reorder(Description, NES),
                                   colour = padj, size = gs_size)) +
                    geom_point() +
                    scale_colour_gradient(low = "#d73027", high = "#4575b4", name = "adj.p") +
                    scale_size_continuous(name = "Set size", range = c(3, 10)) +
                    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
                    labs(title = paste(input$database, "GSEA —", input$contrast),
                         x = "Normalised Enrichment Score", y = NULL) +
                    theme_bw(base_size = 11) +
                    theme(axis.text.y = element_text(size = 9))
            }
        })

        # ---- Ridge Plot -----------------------------------------------------
        output$ridgeplot <- renderPlot({
            obj <- current_gse_obj()
            if (is.null(obj) || !isS4(obj)) {
                return(ggplot() +
                    annotate("text", x = 0.5, y = 0.5,
                             label = "Ridge plot requires enrichResult objects (not available for Hallmarks).\nRun the pipeline to generate RDS files.",
                             hjust = 0.5, vjust = 0.5, size = 4, colour = "grey40") +
                    theme_void())
            }
            tryCatch(
                enrichplot::ridgeplot(obj, showCategory = input$n_terms) +
                    labs(title = paste(input$database, "GSEA —", input$contrast)) +
                    theme_bw(base_size = 11),
                error = function(e) {
                    ggplot() + annotate("text", x=0.5, y=0.5,
                                        label = paste("enrichplot error:", e$message)) + theme_void()
                }
            )
        })

        # ---- Running Sum (gseaplot2) -----------------------------------------
        output$gseaplot <- renderPlot({
            obj <- current_gse_obj()
            if (is.null(obj) || !isS4(obj)) {
                return(ggplot() +
                    annotate("text", x = 0.5, y = 0.5,
                             label = "GSEA running sum plot requires enrichResult objects\n(not available for Hallmarks or when no RDS was saved).",
                             hjust = 0.5, vjust = 0.5, size = 4, colour = "grey40") +
                    theme_void())
            }
            req(input$pathway)
            tryCatch({
                # gseaplot2 takes the index of the gene set in the result object
                df <- as.data.frame(obj)
                idx <- which(df$Description == input$pathway)
                if (length(idx) == 0) idx <- 1L
                enrichplot::gseaplot2(obj, geneSetID = idx,
                                      title = input$pathway)
            }, error = function(e) {
                ggplot() + annotate("text", x=0.5, y=0.5,
                                    label = paste("enrichplot error:", e$message)) + theme_void()
            })
        })

        # ---- Results Table --------------------------------------------------
        output$results_table <- renderDT({
            df <- current_gsea_df()
            if (nrow(df) == 0) return(datatable(data.frame()))
            cols_to_show <- intersect(
                c("Description", "ontology", "NES", "pval", "padj",
                  "size", "setSize", "core_enrichment", "leadingEdge"),
                colnames(df)
            )
            df_show <- df[, cols_to_show, drop = FALSE]
            if ("NES"  %in% colnames(df_show)) df_show$NES  <- round(df_show$NES, 3)
            if ("padj" %in% colnames(df_show)) df_show$padj <- signif(df_show$padj, 3)
            if ("pval" %in% colnames(df_show)) df_show$pval <- signif(df_show$pval, 3)
            datatable(df_show,
                      filter   = "top",
                      options  = list(pageLength = 20, scrollX = TRUE),
                      rownames = FALSE,
                      class    = "compact stripe")
        })
    })
}
