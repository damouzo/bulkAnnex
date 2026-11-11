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
                    numericInput(ns("padj_cut"), "adj.p cutoff",
                                 value = 0.05, min = 0.001, max = 0.5, step = 0.001),
                    numericInput(ns("n_terms"),  "Top N terms",
                                 value = 20, min = 5, max = 50),
                    numericInput(ns("n_gseaplot"), "GSEA plot: top N pathways",
                                 value = 5, min = 1, max = 20),
                    actionButton(ns("apply"), "Apply settings",
                                 icon = icon("play"),
                                 class = "btn-primary w-100 mt-2")
                )
            ),
            column(9,
                # Dynamic tab panel — adds Pathview tab when KEGG is selected
                uiOutput(ns("gsea_panels"))
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

        # ---- Dynamic tab panel ----------------------------------------------
        # Adds Pathview tab only when KEGG database is selected.
        output$gsea_panels <- renderUI({
            db <- input$database %||% "go"
            ns <- session$ns

            base_panels <- list(
                nav_panel("Dot Plot",
                    imageOutput(ns("dotplot"), width = "100%", height = "auto"),
                    div(class = "mt-2",
                        downloadButton(ns("dotplot_dl_png"), "Download PNG",
                                       class = "btn-outline-secondary btn-sm"),
                        downloadButton(ns("dotplot_dl_pdf"), "Download PDF",
                                       class = "btn-outline-secondary btn-sm"))
                ),
                nav_panel("Ridgeplot",
                    imageOutput(ns("ridgeplot"), width = "100%", height = "auto"),
                    div(class = "mt-2",
                        downloadButton(ns("ridgeplot_dl_png"), "Download PNG",
                                       class = "btn-outline-secondary btn-sm"),
                        downloadButton(ns("ridgeplot_dl_pdf"), "Download PDF",
                                       class = "btn-outline-secondary btn-sm"))
                ),
                nav_panel("GSEA Plot",
                    plotOutput(ns("gseaplot"), height = "560px"),
                    div(class = "mt-2",
                        downloadButton(ns("gseaplot_dl_png"), "Download PNG",
                                       class = "btn-outline-secondary btn-sm"),
                        downloadButton(ns("gseaplot_dl_pdf"), "Download PDF",
                                       class = "btn-outline-secondary btn-sm"))
                ),
                nav_panel("Results Table",
                    div(style = "display:flex; flex-direction:column; height:100%;",
                        div(style = "flex:1 1 auto; overflow:auto; min-height:0;",
                            DTOutput(ns("results_table"))),
                        div(style = "flex:0 0 auto; padding-top:8px;",
                            downloadButton(ns("results_dl_csv"), "Download table (CSV)",
                                           class = "btn-outline-secondary btn-sm"))
                    )
                ),
                nav_panel("Common GSEA",
                    fluidRow(
                        column(3,
                            card(
                                card_header("Parameters"),
                                selectInput(ns("cg_db"), "Database",
                                            choices = c("GO" = "go", "KEGG" = "kegg",
                                                        "Reactome" = "reactome"),
                                            selected = "go"),
                                conditionalPanel(
                                    condition = sprintf("input['%s'] == 'go'", ns("cg_db")),
                                    selectInput(ns("cg_go_ont"), "GO Ontology",
                                                choices = c("BP" = "BP", "MF" = "MF", "CC" = "CC"),
                                                selected = "BP")
                                ),
                                numericInput(ns("cg_padj"), "adj.p cutoff",
                                             value = 0.05, min = 0.001, max = 0.5, step = 0.001),
                                radioButtons(ns("cg_type"), "Plot type",
                                             choices = c("Venn Diagram" = "venn",
                                                         "UpSet Plot"   = "upset",
                                                         "Heatmap"      = "heatmap"),
                                             selected = "venn"),
                                conditionalPanel(
                                    condition = sprintf("input['%s'] == 'heatmap'", ns("cg_type")),
                                    selectInput(ns("cg_rowlabels"), "Heatmap row labels",
                                                choices = c("Description (truncated)" = "desc",
                                                            "ID"                      = "id",
                                                            "None"                    = "none"),
                                                selected = "desc")
                                ),
                                hr(),
                                strong("Contrasts to compare"),
                                uiOutput(ns("cg_contrast_check_ui")),
                                hr(),
                                actionButton(ns("cg_apply"), "Apply Changes",
                                             class = "btn-primary w-100", icon = icon("play")),
                                downloadButton(ns("cg_download_png"), "Download PNG",
                                               class = "btn-outline-secondary w-100 mt-2"),
                                downloadButton(ns("cg_download_pdf"), "Download PDF",
                                               class = "btn-outline-secondary w-100 mt-2"),
                                uiOutput(ns("cg_warn"))
                            )
                        ),
                        column(9,
                            card(
                                card_header("Common GSEA"),
                                plotOutput(ns("cg_plot"), height = "620px")
                            )
                        )
                    )
                )
            )

            if (isTRUE(db == "kegg")) {
                base_panels <- c(base_panels, list(
                    nav_panel("Pathview",
                        tags$div(style = "margin-bottom: 0.75rem;",
                            selectInput(ns("pathview_pathway"), "Pathway:",
                                        choices = NULL, width = "100%")
                        ),
                        tags$p(class = "text-muted",
                               tags$em("Top 10 significant KEGG pathways coloured by log2FC.")),
                        imageOutput(ns("pathview_img"), width = "100%", height = "auto"),
                        div(class = "mt-2",
                            downloadButton(ns("pathview_dl_png"), "Download PNG",
                                           class = "btn-outline-secondary btn-sm"))
                    )
                ))
            }

            # TreeDot is always appended as the last tab.
            # Database is taken from the global Settings selector (input$database / input$go_ont).
            base_panels <- c(base_panels, list(
                nav_panel("TreeDot",
                    fluidRow(
                        column(3,
                            card(
                                card_header("TreeDot Settings"),
                                numericInput(ns("td_top_paths"), "Top pathways / contrast",
                                             value = 5, min = 1, max = 20),
                                numericInput(ns("td_clust_num"), "Dendrogram clusters",
                                             value = 3, min = 2, max = 10),
                                tags$hr(),
                                tags$h6("ORA branch labels"),
                                selectInput(ns("td_ora_type"), "ORA database",
                                            choices = c("GO" = "GO", "KEGG" = "KEGG",
                                                        "Reactome" = "Reactome"),
                                            selected = "GO"),
                                conditionalPanel(
                                    condition = sprintf("input['%s'] == 'GO'", ns("td_ora_type")),
                                    selectInput(ns("td_ora_ont"), "ORA GO ontology",
                                                choices = c("BP" = "BP", "MF" = "MF", "CC" = "CC"),
                                                selected = "BP")
                                ),
                                numericInput(ns("td_ora_min_gs"), "ORA min GS size",
                                             value = 10, min = 5, max = 200),
                                numericInput(ns("td_ora_max_gs"), "ORA max GS size",
                                             value = 500, min = 100, max = 2000),
                                numericInput(ns("td_ora_padj"), "ORA adj.p cutoff",
                                             value = 1, min = 0, max = 1, step = 0.05),
                                tags$hr(),
                                strong("Contrasts to include"),
                                uiOutput(ns("td_contrast_ui")),
                                actionButton(ns("td_apply"), "Apply Changes",
                                             icon = icon("play"),
                                             class = "btn-primary w-100 mt-2"),
                                downloadButton(ns("td_download_png"), "Download PNG",
                                               class = "btn-outline-secondary w-100 mt-2"),
                                downloadButton(ns("td_download_pdf"), "Download PDF",
                                               class = "btn-outline-secondary w-100 mt-2")
                            )
                        ),
                        column(9,
                            shinycssloaders::withSpinner(
                                imageOutput(ns("treedot_plot"), width = "100%", height = "auto")
                            )
                        )
                    )
                )
            ))

            do.call(navset_card_tab, base_panels)
        })

        # ---- Frozen display params ------------------------------------------
        plot_params <- reactive({
            list(
                contrast   = input$contrast,
                database   = input$database,
                go_ont     = input$go_ont %||% "BP",
                padj_cut   = input$padj_cut,
                n_terms    = input$n_terms,
                n_gseaplot = input$n_gseaplot
            )
        }) |> bindEvent(input$apply, input$contrast, input$database, input$go_ont,
                        ignoreNULL = TRUE, ignoreInit = FALSE)

        # ---- Resolve the native gseaResult object (or NULL) -----------------
        # gseaResult (S4) available for GO, KEGG and Reactome.
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

        # ---- CSV-based data frame (all tabs except GSEA Plot) ---------------
        current_gsea_df <- reactive({
            req(input$contrast, input$database)
            db_data <- data()$gsea[[input$database]]
            if (is.null(db_data)) return(data.frame())
            df <- db_data[[input$contrast]]
            if (is.null(df) || nrow(df) == 0) return(data.frame())
            # Normalise column names across GO/KEGG/Reactome formats
            if (!"padj"        %in% colnames(df) && "p.adjust"    %in% colnames(df)) df$padj  <- df$p.adjust
            if (!"padj"        %in% colnames(df) && "pval"        %in% colnames(df)) df$padj  <- p.adjust(df$pval, method = "BH")
            if (!"pval"        %in% colnames(df) && "pvalue"      %in% colnames(df)) df$pval  <- df$pvalue
            if (!"Description" %in% colnames(df) && "pathway"     %in% colnames(df)) df$Description <- df$pathway
            if (!"leadingEdge" %in% colnames(df) && "leading_edge"%in% colnames(df)) df$leadingEdge <- df$leading_edge
            df
        })

        # ---- Dot Plot -------------------------------------------------------
        # Priority: S4 object (dynamic, respects n_terms) → CSV fallback → pre-rendered PNG.
        # The pre-rendered PNG is only used as a last resort when no data is in memory.
        # dotplot_state caches the last render (PNG path + ggplot) for the exports.
        dotplot_state <- reactiveVal(NULL)

        output$dotplot <- renderImage({
            p_val <- plot_params()
            obj   <- current_gse_obj()
            df    <- current_gsea_df()
            offline_dbs <- c("kegg", "reactome")
            tmpfile     <- tempfile(fileext = ".png")

            if (!is.null(obj) && isS4(obj)) {
                # Dynamic render — respects n_terms slider
                p <- tryCatch(
                    enrichplot::dotplot(obj, showCategory = p_val$n_terms,
                                        title = paste(toupper(p_val$database),
                                                      "GSEA -", p_val$contrast)) +
                        theme_bw(base_size = 11),
                    error = function(e) {
                        ggplot() +
                            annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                                     label = paste("enrichplot error:", e$message),
                                     hjust = 0.5, vjust = 0.5) +
                            theme_void()
                    }
                )
                ggplot2::ggsave(tmpfile, p, width = 10, height = 7, dpi = 100, bg = "white")
                dotplot_state(list(png = tmpfile, plot = p, w = 10, h = 7))
                return(list(src = tmpfile, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            if (nrow(df) > 0) {
                df_plot <- df %>%
                    filter(padj < p_val$padj_cut) %>%
                    arrange(padj) %>%
                    head(p_val$n_terms) %>%
                    mutate(
                        Description = stringr::str_wrap(Description, 50),
                        gs_size = if ("size"    %in% colnames(.)) size
                                  else if ("setSize" %in% colnames(.)) setSize
                                  else 50
                    )
                if (nrow(df_plot) == 0) {
                    p <- ggplot() +
                        annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                                 label = "No significant terms at current threshold.",
                                 hjust = 0.5, vjust = 0.5) +
                        theme_void()
                } else {
                    p <- ggplot(df_plot, aes(x = NES, y = reorder(Description, NES),
                                             colour = padj, size = gs_size)) +
                        geom_point() +
                        scale_colour_gradient(low = "#d73027", high = "#4575b4", name = "adj.p") +
                        scale_size_continuous(name = "Set size", range = c(3, 10)) +
                        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
                        labs(title = paste(toupper(p_val$database), "GSEA -", p_val$contrast),
                             x = "Normalised Enrichment Score", y = NULL) +
                        theme_bw(base_size = 11) +
                        theme(axis.text.y = element_text(size = 9))
                }
                ggplot2::ggsave(tmpfile, p, width = 10, height = 7, dpi = 100, bg = "white")
                dotplot_state(list(png = tmpfile, plot = p, w = 10, h = 7))
                return(list(src = tmpfile, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            # Last resort: pre-rendered pipeline PNG
            png_path <- get_gsea_png(RESULTS_DIR, p_val$contrast,
                                     p_val$database, p_val$go_ont, type = "dotplot")
            if (!is.null(png_path)) {
                dotplot_state(list(png = png_path, plot = NULL, w = 10, h = 7))
                return(list(src = png_path, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            msg <- if (p_val$database %in% offline_dbs) {
                paste0(toupper(p_val$database), " requires internet access from compute nodes.\n",
                       "Results are unavailable when running on an isolated HPC cluster.")
            } else {
                paste0("No results available for ", toupper(p_val$database), ".")
            }
            p <- ggplot() +
                annotate("text", x = 0.5, y = 0.5, size = 3.5, colour = "grey50",
                         label = msg, hjust = 0.5, vjust = 0.5) +
                theme_void()
            ggplot2::ggsave(tmpfile, p, width = 10, height = 7, dpi = 100, bg = "white")
            dotplot_state(list(png = tmpfile, plot = p, w = 10, h = 7))
            list(src = tmpfile, contentType = "image/png",
                 style = "max-width:100%; height:auto;")
        }, deleteFile = FALSE)

        # ---- Ridgeplot -------------------------------------------------------
        # Priority: S4 object (dynamic) → pre-rendered PNG → error message.
        # IMPORTANT: enrichplot::ridgeplot() can fail on KEGG/Reactome objects
        # when setReadable() has been applied (geneList slot mismatch). In that
        # case we silently fall through to the pre-rendered PNG rather than
        # showing an error — the PNG was generated by the pipeline and is valid.
        # ridgeplot_state caches the last render for the exports.
        ridgeplot_state <- reactiveVal(NULL)

        output$ridgeplot <- renderImage({
            p_val <- plot_params()
            obj   <- current_gse_obj()
            tmpfile <- tempfile(fileext = ".png")

            if (!is.null(obj) && isS4(obj)) {
                df_obj  <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
                col_p   <- if ("p.adjust" %in% colnames(df_obj)) "p.adjust" else "padj"
                n_sig   <- if (col_p %in% colnames(df_obj))
                               sum(!is.na(df_obj[[col_p]]) & df_obj[[col_p]] < p_val$padj_cut)
                           else 0

                ridge_p <- NULL
                if (n_sig == 0) {
                    ridge_p <- ggplot() +
                        annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                                 label = paste0("No significant terms at adj.p < ", p_val$padj_cut),
                                 hjust = 0.5, vjust = 0.5) +
                        theme_void()
                } else {
                    n_show  <- min(p_val$n_terms, n_sig)
                    ridge_p <- tryCatch(
                        enrichplot::ridgeplot(obj, showCategory = n_show) +
                            labs(title = paste(toupper(p_val$database), "Ridgeplot -",
                                               p_val$contrast),
                                 x = "Normalised Enrichment Score") +
                            theme_bw(base_size = 11),
                        error = function(e) NULL   # NULL → fall through to PNG below
                    )
                }

                if (!is.null(ridge_p)) {
                    plot_h <- max(6, min(max(n_sig, 1L) * 0.45 + 2, 16))
                    ggplot2::ggsave(tmpfile, ridge_p, width = 10, height = plot_h,
                                   dpi = 100, bg = "white")
                    ridgeplot_state(list(png = tmpfile, plot = ridge_p, w = 10, h = plot_h))
                    return(list(src = tmpfile, contentType = "image/png",
                                style = "max-width:100%; height:auto;"))
                }
                # S4 ridgeplot failed — fall through to pre-rendered PNG
            }

            # Fallback: pre-rendered pipeline PNG
            png_path <- get_gsea_png(RESULTS_DIR, p_val$contrast,
                                     p_val$database, p_val$go_ont, type = "ridgeplot")
            if (!is.null(png_path)) {
                ridgeplot_state(list(png = png_path, plot = NULL, w = 10, h = 6))
                return(list(src = png_path, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            msg <- if (p_val$database %in% c("kegg", "reactome")) {
                paste0(toupper(p_val$database), " requires internet access \u2014 ridgeplot not available.\n",
                       "Run the pipeline on a node with internet access\n",
                       "to generate the pre-rendered ridgeplot.")
            } else {
                paste0("Ridgeplot not available for ", toupper(p_val$database), ".\n",
                       "Re-run the pipeline to generate it.")
            }
            p <- ggplot() +
                annotate("text", x = 0.5, y = 0.5, size = 3.5, colour = "grey50",
                         label = msg, hjust = 0.5, vjust = 0.5) +
                theme_void()
            ggplot2::ggsave(tmpfile, p, width = 10, height = 6, dpi = 100, bg = "white")
            ridgeplot_state(list(png = tmpfile, plot = p, w = 10, h = 6))
            list(src = tmpfile, contentType = "image/png",
                 style = "max-width:100%; height:auto;")
        }, deleteFile = FALSE)

        # Dot Plot / Ridgeplot exports: reuse the cached render (ggplot when
        # available, otherwise the pre-rendered PNG for PNG downloads).
        dl_from_state <- function(state, file, ext) {
            if (is.null(state)) req(FALSE)
            if (!is.null(state$plot)) {
                ggplot2::ggsave(file, state$plot, device = ext,
                                width = state$w %||% 10, height = state$h %||% 7,
                                dpi = 300, bg = "white")
            } else if (ext == "png" && !is.null(state$png) && file.exists(state$png)) {
                file.copy(state$png, file, overwrite = TRUE)
            } else {
                req(FALSE)
            }
        }

        output$dotplot_dl_png <- downloadHandler(
            filename = function() paste0(input$contrast %||% "contrast", "_dotplot.png"),
            content  = function(file) dl_from_state(dotplot_state(), file, "png")
        )

        output$dotplot_dl_pdf <- downloadHandler(
            filename = function() paste0(input$contrast %||% "contrast", "_dotplot.pdf"),
            content  = function(file) dl_from_state(dotplot_state(), file, "pdf")
        )

        output$ridgeplot_dl_png <- downloadHandler(
            filename = function() paste0(input$contrast %||% "contrast", "_ridgeplot.png"),
            content  = function(file) dl_from_state(ridgeplot_state(), file, "png")
        )

        output$ridgeplot_dl_pdf <- downloadHandler(
            filename = function() paste0(input$contrast %||% "contrast", "_ridgeplot.pdf"),
            content  = function(file) dl_from_state(ridgeplot_state(), file, "pdf")
        )

        # ---- GSEA Plot (gseaplot2 — top N pathways by NES) ------------------
        # Requires a gseaResult (S4) object: available for GO, KEGG and Reactome.
        # Build the gseaplot2 object; shared by renderPlot and the exports.
        gseaplot_obj <- function() {
            p_val <- plot_params()
            obj   <- current_gse_obj()

            if (is.null(obj) || !isS4(obj)) {
                msg <- if (p_val$database %in% c("kegg", "reactome")) {
                    paste0(toupper(p_val$database), " results unavailable (requires internet).\n",
                           "Run the pipeline on a node with internet access\n",
                           "to enable this plot.")
                } else {
                    "No gseaResult object available for this contrast."
                }
                return(ggplot() +
                    annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                             label = msg, hjust = 0.5, vjust = 0.5) +
                    theme_void())
            }

            df_obj <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
            if (nrow(df_obj) == 0) {
                return(ggplot() +
                    annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                             label = "No significant gene sets found.",
                             hjust = 0.5, vjust = 0.5) +
                    theme_void())
            }

            n_top <- min(p_val$n_gseaplot, nrow(df_obj))
            tryCatch(
                enrichplot::gseaplot2(obj, geneSetID = seq_len(n_top),
                                      title = paste0(toupper(p_val$database),
                                                     " - top ", n_top, " pathways")),
                error = function(e) {
                    ggplot() +
                        annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                                 label = paste("gseaplot2 error:", e$message),
                                 hjust = 0.5, vjust = 0.5) +
                        theme_void()
                }
            )
        }

        output$gseaplot <- renderPlot({ gseaplot_obj() })

        output$gseaplot_dl_png <- downloadHandler(
            filename = function() paste0(input$contrast %||% "contrast", "_gseaplot.png"),
            content  = function(file) {
                ggplot2::ggsave(file, gseaplot_obj(), device = "png",
                                width = 10, height = 7, dpi = 300, bg = "white")
            }
        )

        output$gseaplot_dl_pdf <- downloadHandler(
            filename = function() paste0(input$contrast %||% "contrast", "_gseaplot.pdf"),
            content  = function(file) {
                ggplot2::ggsave(file, gseaplot_obj(), device = "pdf",
                                width = 10, height = 7, bg = "white")
            }
        )

        # ---- Results Table --------------------------------------------------
        gsea_results_table_df <- reactive({
            p_val <- plot_params()
            df    <- current_gsea_df()
            if (nrow(df) == 0) return(data.frame())
            if (p_val$database == "go" && "ontology" %in% colnames(df)) {
                df <- df[df$ontology == p_val$go_ont, , drop = FALSE]
            }
            cols_to_show <- intersect(
                c("Description", "ontology", "NES", "pval", "padj",
                  "size", "setSize", "core_enrichment", "leadingEdge"),
                colnames(df)
            )
            df_show <- df[, cols_to_show, drop = FALSE]
            if ("NES"  %in% colnames(df_show)) df_show$NES  <- round(df_show$NES,  3)
            if ("padj" %in% colnames(df_show)) df_show$padj <- signif(df_show$padj, 3)
            if ("pval" %in% colnames(df_show)) df_show$pval <- signif(df_show$pval, 3)
            df_show
        })

        output$results_table <- renderDT({
            df_show <- gsea_results_table_df()
            if (nrow(df_show) == 0) return(datatable(data.frame()))

            le_col <- intersect(c("core_enrichment", "leadingEdge"), colnames(df_show))
            col_defs <- if (length(le_col) > 0) {
                le_idx <- which(colnames(df_show) == le_col[1]) - 1L
                list(list(
                    targets = le_idx,
                    render  = DT::JS(
                        "function(data, type, row) {
                            if (type !== 'display' || !data || data.length <= 60)
                                return data || '';
                            return '<details><summary style=\"cursor:pointer;color:#2c7bb6\">' +
                                   data.substring(0, 60) + '\u2026</summary>' +
                                   '<span style=\"word-break:break-all;font-size:0.85em\">' +
                                   data + '</span></details>';
                        }"
                    )
                ))
            } else list()

            datatable(df_show,
                      filter      = "top",
                      options     = list(pageLength = 20, scrollX = TRUE,
                                         columnDefs = col_defs),
                      rownames    = FALSE,
                      escape      = FALSE,
                      class       = "compact stripe")
        })

        output$results_dl_csv <- downloadHandler(
            filename = function() {
                p_val <- plot_params()
                paste0(p_val$contrast %||% "contrast", "_",
                       toupper(p_val$database %||% ""), "_GSEA_results.csv")
            },
            content = function(file) {
                df_show <- gsea_results_table_df()
                req(nrow(df_show) > 0)
                write.csv(df_show, file, row.names = FALSE)
            }
        )

        # =====================================================================
        # COMMON GSEA — cross-contrast pathway comparison (Venn / UpSet / Heatmap)
        # Compares significant pathways (by ID, value = NES) across contrasts,
        # independently of the global Settings selectors above.
        # =====================================================================
        cg_pkgs_ok <- all(sapply(
            c("ggVennDiagram", "ComplexUpset", "ComplexHeatmap", "circlize"),
            requireNamespace, quietly = TRUE
        ))

        # ---- Populate contrast checkboxes -----------------------------------
        output$cg_contrast_check_ui <- renderUI({
            ns <- session$ns
            cids <- unique(unlist(lapply(data()$gsea, names)))
            if (length(cids) == 0) cids <- names(data()$dge)
            if (length(cids) == 0) {
                return(p("No GSEA results available.", style = "color:grey;font-size:0.9em;"))
            }
            checkboxGroupInput(ns("cg_contrasts"), label = NULL, choices = cids, selected = cids)
        })

        output$cg_warn <- renderUI({
            req(input$cg_type == "venn")
            n <- length(input$cg_contrasts)
            if (!is.null(n) && n > 5) {
                div(class = "alert alert-warning mt-2", style = "font-size:0.83em; padding:8px;",
                    tags$strong("\u26a0 Too many contrasts for a Venn Diagram (max 5)."),
                    tags$br(),
                    "Please select \u22645 contrasts or switch to UpSet Plot.")
            }
        })

        # Normalise a GSEA CSV data.frame to canonical columns.
        cg_normalise_df <- function(df, db, ont) {
            if (is.null(df) || nrow(df) == 0) return(NULL)
            if (!"padj" %in% colnames(df) && "p.adjust" %in% colnames(df)) df$padj <- df$p.adjust
            if (!"padj" %in% colnames(df) && "pval"     %in% colnames(df)) df$padj <- p.adjust(df$pval, method = "BH")
            if (!"ID"   %in% colnames(df)) df$ID <- df$Description
            if (!"NES"  %in% colnames(df)) df$NES <- NA_real_
            if (db == "go" && "ontology" %in% colnames(df)) {
                df <- df[toupper(df$ontology) == ont, , drop = FALSE]
            }
            df
        }

        cg_sets <- eventReactive(input$cg_apply, {
            db       <- input$cg_db %||% "go"
            ont      <- toupper(input$cg_go_ont %||% "BP")
            padj_thr <- input$cg_padj %||% 0.05
            selected <- input$cg_contrasts
            validate(need(length(selected) >= 2,
                          "Please select at least 2 contrasts to compare."))

            base <- data()$gsea[[db]]
            if (is.null(base)) base <- list()

            dfs <- setNames(lapply(selected, function(cid) {
                df <- cg_normalise_df(base[[cid]], db, ont)
                if (is.null(df)) return(NULL)
                df[!is.na(df$padj) & df$padj < padj_thr, , drop = FALSE]
            }), selected)

            sets <- lapply(dfs, function(d) if (is.null(d) || nrow(d) == 0) character(0) else unique(d$ID))

            all_df <- do.call(rbind, lapply(dfs, function(d) {
                if (is.null(d) || nrow(d) == 0) return(NULL)
                d[, c("ID", "Description"), drop = FALSE]
            }))
            labels <- if (!is.null(all_df) && nrow(all_df) > 0) {
                all_df <- all_df[!duplicated(all_df$ID), , drop = FALSE]
                setNames(all_df$Description, all_df$ID)
            } else character(0)

            list(sets = sets, labels = labels, dfs = dfs, db = db, ont = ont)
        }, ignoreNULL = FALSE)

        cg_trunc <- function(x, n = 45) {
            ifelse(is.na(x), "", ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x))
        }

        cg_draw_venn <- function(sets) {
            if (length(sets) > 5) {
                grid::grid.text("Too many contrasts for a Venn Diagram (max 5).",
                                gp = grid::gpar(fontsize = 12, col = "firebrick"))
                return(invisible(NULL))
            }
            if (all(lengths(sets) == 0)) {
                grid::grid.text("No significant pathways found.",
                                gp = grid::gpar(fontsize = 12, col = "grey50"))
                return(invisible(NULL))
            }
            p <- ggVennDiagram::ggVennDiagram(
                    sets, label = "both", set_color = rep("black", length(sets))) +
                ggplot2::scale_fill_gradient(low = "white", high = "#4d9221") +
                ggplot2::scale_color_manual(values = rep("black", length(sets))) +
                ggplot2::ggtitle("Common GSEA pathways") +
                ggplot2::theme(
                    legend.position = "none",
                    text = ggplot2::element_text(color = "black"),
                    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
                    plot.margin = ggplot2::margin(16, 36, 16, 36, unit = "pt")) +
                ggplot2::coord_cartesian(clip = "off")
            print(p)
        }

        cg_draw_upset <- function(sets) {
            all_paths <- unique(unlist(sets))
            if (length(all_paths) == 0) {
                grid::grid.text("No significant pathways found.",
                                gp = grid::gpar(fontsize = 12, col = "grey50"))
                return(invisible(NULL))
            }
            mat <- as.data.frame(
                sapply(names(sets), function(cid) all_paths %in% sets[[cid]]),
                stringsAsFactors = FALSE)
            mat[] <- lapply(mat, as.logical)
            p <- ComplexUpset::upset(
                    mat, intersect = names(sets), name = "Contrast",
                    set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(fill = "#4d9221")),
                    width_ratio = 0.2, sort_intersections_by = "degree") +
                ggplot2::ggtitle("Common GSEA pathways") +
                ggplot2::theme(
                    text = ggplot2::element_text(color = "black"),
                    axis.text = ggplot2::element_text(color = "black"),
                    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13))
            print(p)
        }

        cg_draw_heatmap <- function(sets, labels, dfs, labels_mode) {
            common <- Reduce(intersect, sets)
            if (length(common) < 2) {
                grid::grid.text(
                    paste0("No pathways significant in ALL selected contrasts (n = ",
                           length(common), ")."),
                    gp = grid::gpar(fontsize = 12, col = "grey50"))
                return(invisible(NULL))
            }
            nes_mat <- sapply(names(dfs), function(cid) {
                d <- dfs[[cid]]
                if (is.null(d) || nrow(d) == 0) return(rep(NA_real_, length(common)))
                setNames(d$NES, d$ID)[common]
            })
            nes_mat <- as.matrix(nes_mat)
            rownames(nes_mat) <- common

            show_rn <- nrow(nes_mat) <= 50
            if (labels_mode == "none") {
                show_rn <- FALSE
            } else if (labels_mode == "id") {
                rownames(nes_mat) <- make.unique(common)
            } else {
                labs <- labels[common]
                miss <- is.na(labs) | !nzchar(labs)
                labs[miss] <- common[miss]
                rownames(nes_mat) <- make.unique(cg_trunc(labs, 45))
            }

            lim     <- max(stats::quantile(abs(nes_mat), 0.98, na.rm = TRUE), 1e-6)
            col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))
            ht <- ComplexHeatmap::Heatmap(
                nes_mat,
                name               = "NES",
                col                = col_fun,
                na_col             = "#cccccc",
                show_row_names     = show_rn,
                row_names_gp       = grid::gpar(fontsize = 7),
                row_names_max_width= grid::unit(7, "cm"),
                column_names_gp    = grid::gpar(fontsize = 9),
                cluster_rows       = TRUE,
                cluster_columns    = FALSE,
                column_title       = paste0("Common GSEA pathways (n = ", nrow(nes_mat),
                                            " x ", ncol(nes_mat), " contrasts)"),
                column_title_gp    = grid::gpar(fontsize = 12, fontface = "bold"),
                heatmap_legend_param = list(title_gp  = grid::gpar(fontsize = 9),
                                            labels_gp = grid::gpar(fontsize = 8))
            )
            ComplexHeatmap::draw(ht)
        }

        cg_render_plot <- function() {
            if (!cg_pkgs_ok) {
                grid::grid.text("Package(s) unavailable. Please rebuild the container.",
                                gp = grid::gpar(fontsize = 12, col = "firebrick"))
                return(invisible(NULL))
            }
            d <- cg_sets()
            req(d)
            if (length(d$sets) < 2) {
                grid::grid.text("Select at least 2 contrasts.",
                                gp = grid::gpar(fontsize = 13, col = "grey50"))
                return(invisible(NULL))
            }
            type <- input$cg_type %||% "venn"
            if (type == "venn") {
                cg_draw_venn(d$sets)
            } else if (type == "upset") {
                cg_draw_upset(d$sets)
            } else {
                cg_draw_heatmap(d$sets, d$labels, d$dfs, input$cg_rowlabels %||% "desc")
            }
        }

        output$cg_plot <- renderPlot({ cg_render_plot() })

        output$cg_download_png <- downloadHandler(
            filename = function() paste0("common_gsea_", input$cg_type %||% "plot", ".png"),
            content  = function(file) {
                png(file, width = 1600, height = 1200, res = 200)
                on.exit(dev.off(), add = TRUE)
                cg_render_plot()
            }
        )

        output$cg_download_pdf <- downloadHandler(
            filename = function() paste0("common_gsea_", input$cg_type %||% "plot", ".pdf"),
            content  = function(file) {
                pdf(file, width = 10, height = 8)
                on.exit(dev.off(), add = TRUE)
                cg_render_plot()
            }
        )

        # ---- Pathview (KEGG only) -------------------------------------------
        # Populate pathway dropdown when KEGG contrast changes.
        observeEvent(list(input$contrast, input$database), {
            req(input$database == "kegg", input$contrast)
            pv_files <- get_pathview_files(RESULTS_DIR, input$contrast)
            if (length(pv_files) == 0) {
                updateSelectInput(session, "pathview_pathway",
                                  choices = c("No pathview files available" = ""))
                return()
            }
            # Enrich labels with pathway names from KEGG CSV if available
            df_kegg <- isolate(current_gsea_df())
            pw_ids  <- names(pv_files)
            if (nrow(df_kegg) > 0 && all(c("ID", "Description") %in% colnames(df_kegg))) {
                name_map <- setNames(df_kegg$Description, df_kegg$ID)
                labels   <- ifelse(pw_ids %in% names(name_map),
                                   paste0(pw_ids, " — ", name_map[pw_ids]),
                                   pw_ids)
            } else {
                labels <- pw_ids
            }
            choices <- setNames(pw_ids, labels)
            updateSelectInput(session, "pathview_pathway",
                              choices = choices, selected = choices[1])
        }, ignoreNULL = TRUE)

        # Display the selected pathview image.
        output$pathview_img <- renderImage({
            req(input$database == "kegg")
            pw_id <- input$pathview_pathway
            if (is.null(pw_id) || nchar(pw_id) == 0) {
                tmpfile <- tempfile(fileext = ".png")
                p <- ggplot() +
                    annotate("text", x = 0.5, y = 0.5, size = 3.5, colour = "grey50",
                             label = paste0(
                                 "No Pathview images found for this contrast.\n\n",
                                 "Pathview downloads KEGG pathway maps at runtime\n",
                                 "and requires full internet access from compute nodes.\n\n",
                                 "KEGG GSEA ran successfully (KEGG API was reachable),\n",
                                 "but Pathview uses a different download mechanism\n",
                                 "that may be blocked by the HPC firewall.\n\n",
                                 "To generate Pathview images, re-run the pipeline on a\n",
                                 "node with unrestricted internet access, or contact\n",
                                 "your HPC team to allow outbound HTTP to kegg.jp."
                             ),
                             hjust = 0.5, vjust = 0.5) +
                    theme_void()
                ggplot2::ggsave(tmpfile, p, width = 8, height = 5, dpi = 100, bg = "white")
                return(list(src = tmpfile, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }
            pv_files <- get_pathview_files(RESULTS_DIR, input$contrast)
            f <- pv_files[[pw_id]]
            if (is.null(f) || !file.exists(f)) {
                tmpfile <- tempfile(fileext = ".png")
                p <- ggplot() +
                    annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                             label = paste0("File not found for pathway: ", pw_id),
                             hjust = 0.5, vjust = 0.5) +
                    theme_void()
                ggplot2::ggsave(tmpfile, p, width = 8, height = 5, dpi = 100, bg = "white")
                return(list(src = tmpfile, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }
            list(src = f, contentType = "image/png",
                 style = "max-width:100%; height:auto;")
        }, deleteFile = FALSE)

        output$pathview_dl_png <- downloadHandler(
            filename = function() {
                pw <- input$pathview_pathway %||% "pathview"
                paste0(input$contrast %||% "contrast", "_pathview_", pw, ".png")
            },
            content = function(file) {
                req(input$database == "kegg")
                pw <- input$pathview_pathway
                req(!is.null(pw), nchar(pw) > 0)
                f <- get_pathview_files(RESULTS_DIR, input$contrast)[[pw]]
                req(!is.null(f), file.exists(f))
                file.copy(f, file, overwrite = TRUE)
            }
        )

        # ---- TreeDot --------------------------------------------------------
        # td_gsea_db tag is derived from the global Settings selectors:
        #   input$database ("go"/"kegg"/"reactome") + input$go_ont ("BP"/"MF"/"CC")
        # The user may also restrict which contrasts are included in the plot.
        td_applied <- reactiveVal(NULL)   # list(tag, contrasts) frozen on Apply
        td_cache   <- reactiveVal(NULL)   # list(plot, w, h, tag) of the last render

        # Contrasts available for the TreeDot (from DGE results).
        output$td_contrast_ui <- renderUI({
            ns   <- session$ns
            cids <- names(data()$dge)
            if (length(cids) == 0) {
                return(p("No DGE results available.", style = "color:grey;font-size:0.9em;"))
            }
            checkboxGroupInput(ns("td_contrasts"), label = NULL,
                               choices = cids, selected = cids)
        })

        td_selected_contrasts <- reactive({
            all_c <- names(data()$dge)
            sel   <- input$td_contrasts
            if (is.null(sel)) all_c else intersect(all_c, sel)
        })

        # Reset the applied state whenever database or ontology changes so the
        # pre-rendered PNG is shown first; computation only fires on Apply.
        observeEvent(list(input$database, input$go_ont), {
            td_applied(NULL); td_cache(NULL)
        }, ignoreInit = TRUE)

        observeEvent(input$td_apply, {
            db  <- input$database %||% "go"
            ont <- input$go_ont   %||% "BP"
            tag <- switch(db,
                "go"       = paste0("GO_", toupper(ont)),
                "kegg"     = "KEGG",
                "reactome" = "Reactome",
                "GO_BP"
            )
            td_applied(list(tag = tag, contrasts = sort(td_selected_contrasts())))
            td_cache(NULL)
        })

        # Small placeholder PNG for the various "nothing to show" states.
        td_placeholder <- function(msg, w = 10, h = 6) {
            tmpfile <- tempfile(fileext = ".png")
            p <- ggplot() +
                annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                         size = 3.7, colour = "grey50", label = msg) +
                theme_void()
            ggplot2::ggsave(tmpfile, p, width = w, height = h, dpi = 100, bg = "white")
            tmpfile
        }

        output$treedot_plot <- renderImage({
            db  <- input$database %||% "go"
            ont <- input$go_ont   %||% "BP"
            td_gsea_db <- switch(db,
                "go"       = paste0("GO_", toupper(ont)),
                "kegg"     = "KEGG",
                "reactome" = "Reactome",
                "GO_BP"
            )
            gsea_db  <- switch(db, "go" = "GO", "kegg" = "KEGG", "Reactome")
            gsea_ont <- if (db == "go") toupper(ont) else NA_character_

            all_contrasts <- names(data()$dge)
            sel           <- sort(td_selected_contrasts())
            all_selected  <- length(all_contrasts) > 0 && setequal(sel, all_contrasts)

            applied      <- td_applied()
            apply_active <- !is.null(applied) && identical(applied$tag, td_gsea_db) &&
                            identical(applied$contrasts, sel)

            # 1. No Apply yet → show the pre-rendered PNG only when every contrast
            #    is selected, since that is what the pre-rendered PNG contains.
            if (!apply_active) {
                if (all_selected && length(sel) >= 2) {
                    png_path <- get_treedot_png(RESULTS_DIR, td_gsea_db)
                    if (!is.null(png_path)) {
                        return(list(src = png_path, contentType = "image/png",
                                    style = "max-width:100%; height:auto;"))
                    }
                }
                msg <- if (length(sel) < 2)
                           "Select at least 2 contrasts for the TreeDot."
                       else if (!all_selected)
                           paste0("Contrast subset selected.\n",
                                  "Click 'Apply Changes' to compute the TreeDot.")
                       else paste0("No pre-rendered TreeDot found for ", td_gsea_db, ".\n\n",
                                   "Run the pipeline to generate it, or\n",
                                   "click 'Apply Changes' to compute now.")
                return(list(src = td_placeholder(msg), contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            if (length(sel) < 2) {
                return(list(src = td_placeholder(
                            "Need >= 2 contrasts in DGE results\nto compute the TreeDot."),
                            contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            # Detect organism from the first available RDS
            rds_list  <- data()$treedot_rds
            first_cmp <- if (!is.null(rds_list) && length(rds_list) > 0) rds_list[[1]] else NULL
            call_str  <- tryCatch(paste(deparse(first_cmp@.call), collapse = " "),
                                  error = function(e) "")
            org_db_str <- if (grepl("org\\.Mm", call_str)) "org.Mm.eg.db" else "org.Hs.eg.db"
            kegg_org   <- if (org_db_str == "org.Mm.eg.db") "mmu" else "hsa"

            # Reuse the stored RDS only for the full contrast set; otherwise
            # recompute compareCluster for the selected subset (~minutes).
            use_rds <- all_selected && !is.null(rds_list) && !is.null(rds_list[[td_gsea_db]])
            cmp <- if (use_rds) {
                rds_list[[td_gsea_db]]
            } else {
                dge_list <- data()$dge[sel]
                tryCatch(
                    build_treedot_cmp(dge_list, org_db_str, kegg_org,
                                      gsea_db, gsea_ont %||% "BP",
                                      min_gs = 15L, max_gs = 500L),
                    error = function(e) { message("build_treedot_cmp: ", e$message); NULL }
                )
            }

            if (is.null(cmp)) {
                return(list(src = td_placeholder(paste0(
                            "compareCluster failed.\n\n",
                            "Check GSEA results are available\n",
                            "and internet is accessible (KEGG).")),
                            contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            # Visual params — frozen via isolate so they only change on Apply.
            top_paths  <- isolate(input$td_top_paths)  %||% 5L
            clust_num  <- isolate(input$td_clust_num)  %||% 3L
            ora_type   <- isolate(input$td_ora_type)   %||% "GO"
            ora_ont    <- isolate(input$td_ora_ont)    %||% "BP"
            ora_min_gs <- isolate(input$td_ora_min_gs) %||% 10L
            ora_max_gs <- isolate(input$td_ora_max_gs) %||% 500L
            ora_padj   <- isolate(input$td_ora_padj)   %||% 1.0

            cmp_info <- get_rds_db_info(cmp)
            keytype  <- cmp_info$keytype

            fort_tmp    <- tryCatch(
                fortify(cmp, showCategory = top_paths, includeAll = TRUE, split = NULL),
                error = function(e) data.frame(Description = character(0), Cluster = character(0))
            )
            n_paths     <- max(length(unique(fort_tmp$Description)), 1L)
            n_contrasts <- max(length(unique(sub("\n.*", "", fort_tmp$Cluster))), 1L)
            plot_h      <- max(8, n_paths * 0.35 + 4)
            plot_w      <- max(12, n_contrasts * 1.2 + 6)

            p <- tryCatch(
                plot_treedot(cmp        = cmp,
                             top_paths  = top_paths,
                             clust_num  = clust_num,
                             ora_type   = ora_type,
                             ora_ont    = ora_ont,
                             ora_min_gs = ora_min_gs,
                             ora_max_gs = ora_max_gs,
                             ora_padj   = ora_padj,
                             org_db_str = org_db_str,
                             kegg_org   = kegg_org,
                             keytype    = keytype),
                error = function(e) {
                    ggplot() +
                        annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                                 size = 4, colour = "grey40",
                                 label = paste("TreeDot error:", e$message)) +
                        theme_void()
                }
            )
            td_cache(list(plot = p, w = plot_w, h = plot_h, tag = td_gsea_db))

            tmpfile <- tempfile(fileext = ".png")
            ggplot2::ggsave(tmpfile, p, width = plot_w, height = plot_h, dpi = 100, bg = "white")
            list(src = tmpfile, contentType = "image/png",
                 style = "max-width:100%; height:auto;")
        }, deleteFile = FALSE)

        # Shared TreeDot export: reuse the last rendered plot when available,
        # otherwise fall back to the pre-rendered pipeline file.
        td_export <- function(file, ext) {
            c <- td_cache()
            if (!is.null(c) && !is.null(c$plot)) {
                ggplot2::ggsave(file, c$plot, device = ext,
                                width = c$w, height = c$h, dpi = 300, bg = "white")
                return(invisible(NULL))
            }
            db  <- input$database %||% "go"
            ont <- input$go_ont %||% "BP"
            tag <- switch(db,
                "go"       = paste0("GO_", toupper(ont)),
                "kegg"     = "KEGG",
                "reactome" = "Reactome",
                "GO_BP"
            )
            cand <- file.path(RESULTS_DIR, "gsea", "treedot",
                              paste0("gsea_treedot_", tag, ".", ext))
            if (file.exists(cand)) {
                file.copy(cand, file, overwrite = TRUE)
            } else {
                req(FALSE)
            }
        }

        output$td_download_png <- downloadHandler(
            filename = function() paste0("gsea_treedot_", input$database %||% "go", ".png"),
            content  = function(file) td_export(file, "png")
        )

        output$td_download_pdf <- downloadHandler(
            filename = function() paste0("gsea_treedot_", input$database %||% "go", ".pdf"),
            content  = function(file) td_export(file, "pdf")
        )
    })
}
