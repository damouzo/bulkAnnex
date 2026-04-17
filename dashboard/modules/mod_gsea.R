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
                    imageOutput(ns("dotplot"), width = "100%", height = "auto")
                ),
                nav_panel("Ridgeplot",
                    imageOutput(ns("ridgeplot"), width = "100%", height = "auto")
                ),
                nav_panel("GSEA Plot",
                    plotOutput(ns("gseaplot"), height = "560px")
                ),
                nav_panel("Results Table",
                    DTOutput(ns("results_table"))
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
                        imageOutput(ns("pathview_img"), width = "100%", height = "auto")
                    )
                ))
            }

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
        # Serves the pre-rendered pipeline PNG; falls back to dynamic ggplot.
        output$dotplot <- renderImage({
            p_val <- plot_params()

            png_path <- get_gsea_png(RESULTS_DIR, p_val$contrast,
                                     p_val$database, p_val$go_ont, type = "dotplot")
            if (!is.null(png_path)) {
                return(list(src = png_path, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            # Dynamic fallback
            obj <- current_gse_obj()
            df  <- current_gsea_df()
            offline_dbs <- c("kegg", "reactome")
            tmpfile     <- tempfile(fileext = ".png")

            if (!is.null(obj) && isS4(obj)) {
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
            } else if (nrow(df) > 0) {
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
            } else {
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
            }

            ggplot2::ggsave(tmpfile, p, width = 10, height = 7, dpi = 100, bg = "white")
            list(src = tmpfile, contentType = "image/png",
                 style = "max-width:100%; height:auto;")
        }, deleteFile = FALSE)

        # ---- Ridgeplot -------------------------------------------------------
        # Shows pre-rendered pipeline ridgeplot PNG.
        # Fallback: enrichplot::ridgeplot() for S4 objects (GO/KEGG/Reactome).
        output$ridgeplot <- renderImage({
            p_val <- plot_params()

            png_path <- get_gsea_png(RESULTS_DIR, p_val$contrast,
                                     p_val$database, p_val$go_ont, type = "ridgeplot")
            if (!is.null(png_path)) {
                return(list(src = png_path, contentType = "image/png",
                            style = "max-width:100%; height:auto;"))
            }

            # Dynamic fallback for S4-based databases (GO/KEGG/Reactome)
            tmpfile <- tempfile(fileext = ".png")
            obj     <- current_gse_obj()

            if (!is.null(obj) && isS4(obj)) {
                df_obj  <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
                col_p   <- if ("p.adjust" %in% colnames(df_obj)) "p.adjust" else "padj"
                n_sig   <- if (col_p %in% colnames(df_obj))
                               sum(!is.na(df_obj[[col_p]]) & df_obj[[col_p]] < p_val$padj_cut)
                           else 0
                if (n_sig == 0) {
                    p <- ggplot() +
                        annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                                 label = paste0("No significant terms at adj.p < ", p_val$padj_cut),
                                 hjust = 0.5, vjust = 0.5) +
                        theme_void()
                } else {
                    n_show <- min(p_val$n_terms, n_sig)
                    p <- tryCatch(
                        enrichplot::ridgeplot(obj, showCategory = n_show) +
                            labs(title = paste(toupper(p_val$database), "Ridgeplot -",
                                               p_val$contrast)) +
                            theme_bw(base_size = 11),
                        error = function(e) {
                            ggplot() +
                                annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey40",
                                         label = paste("ridgeplot error:", e$message),
                                         hjust = 0.5, vjust = 0.5) +
                                theme_void()
                        }
                    )
                }
                plot_h <- max(6, min(n_sig * 0.45 + 2, 16))
            } else {
                msg <- if (p_val$database %in% c("kegg", "reactome")) {
                    paste0(toupper(p_val$database), " requires internet access — ridgeplot not available.\n",
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
                plot_h <- 6
            }

            ggplot2::ggsave(tmpfile, p, width = 10, height = plot_h, dpi = 100, bg = "white")
            list(src = tmpfile, contentType = "image/png",
                 style = "max-width:100%; height:auto;")
        }, deleteFile = FALSE)

        # ---- GSEA Plot (gseaplot2 — top N pathways by NES) ------------------
        # Requires a gseaResult (S4) object: available for GO, KEGG and Reactome.
        output$gseaplot <- renderPlot({
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
        })

        # ---- Results Table --------------------------------------------------
        output$results_table <- renderDT({
            p_val <- plot_params()
            df    <- current_gsea_df()
            if (nrow(df) == 0) return(datatable(data.frame()))
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
                    annotate("text", x = 0.5, y = 0.5, size = 4, colour = "grey50",
                             label = paste0("No Pathview images found.\n\n",
                                            "Pathview requires internet access from compute nodes.\n",
                                            "Run the pipeline on a node with KEGG API access."),
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
    })
}
