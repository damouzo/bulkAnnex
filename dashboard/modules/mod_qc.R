# modules/mod_qc.R — QC tab

mod_qc_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12, h3("Quality Control"))
        ),
        fluidRow(
            column(4,
                card(
                    card_header("QC Metrics"),
                    DTOutput(ns("metrics_table"))
                )
            ),
            column(8,
                navset_card_tab(
                    nav_panel("Library Sizes",
                        plotlyOutput(ns("plot_library_sizes"), height = "auto")
                    ),
                    nav_panel("Count Distribution",
                        imageOutput(ns("plot_count_dist"), width = "100%", height = "auto")
                    ),
                    nav_panel("PCA",
                        fluidRow(
                            column(3,
                                tags$div(
                                    style = "border-right:1px solid #dee2e6; min-height:480px; padding:6px 14px 6px 4px;",
                                    tags$p(tags$strong("PCA options"),
                                           style = "font-size:0.9rem; margin-bottom:10px;"),
                                    selectInput(ns("pca_color_by"), "Color by:",
                                        choices  = c("Condition" = "condition"),
                                        selected = "condition",
                                        width    = "100%"
                                    ),
                                    tags$hr(style = "margin:8px 0;"),
                                    tags$p("Samples to include:",
                                           style = "font-size:0.83rem; margin-bottom:3px;"),
                                    selectizeInput(ns("pca_samples"), label = NULL,
                                        choices  = c(),
                                        multiple = TRUE,
                                        width    = "100%",
                                        options  = list(
                                            placeholder = "All samples",
                                            plugins     = list("remove_button")
                                        )
                                    ),
                                    actionButton(ns("pca_apply"), "Apply",
                                        class = "btn-primary btn-sm",
                                        style = "width:100%; margin-top:8px;"),
                                    actionButton(ns("pca_reset"), "Reset to all",
                                        class = "btn-outline-secondary btn-sm",
                                        style = "width:100%; margin-top:4px;")
                                )
                            ),
                            column(9,
                                plotlyOutput(ns("plot_pca"), height = "480px")
                            )
                        )
                    ),
                    nav_panel("Correlation Heatmap",
                        imageOutput(ns("corr_heatmap"), width = "100%", height = "auto")
                    )
                )
            )
        )
    )
}

mod_qc_server <- function(id, app_data) {
    moduleServer(id, function(input, output, session) {
        data <- reactive(app_data())

        qc_metrics <- reactive({
            data()$qc_metrics
        })

        output$metrics_table <- renderDT({
            df <- qc_metrics()
            if (nrow(df) == 0) return(datatable(data.frame(Message = "No QC metrics available")))
            if ("library_size" %in% colnames(df))
                df$library_size <- format(df$library_size, big.mark = ",", scientific = FALSE)
            if ("n_detected" %in% colnames(df))
                df$n_detected   <- format(df$n_detected,   big.mark = ",", scientific = FALSE)
            datatable(df,
                      options  = list(pageLength = 25,
                                      scrollX    = TRUE,
                                      scrollY    = "450px",
                                      dom        = "frtip"),
                      rownames = FALSE,
                      class    = "compact stripe hover")
        })

        output$plot_library_sizes <- renderPlotly({
            df <- qc_metrics()
            if (nrow(df) == 0) return(plot_ly() %>% layout(title = "No data"))
            df$sample <- factor(df$sample, levels = df$sample[order(df$library_size)])

            # Generate as many perceptually distinct colours as there are conditions
            conditions  <- unique(df$condition)
            n_cond      <- length(conditions)
            pal         <- setNames(
                hcl.colors(n_cond, palette = "Dark 3", fixup = TRUE),
                conditions
            )

            # Scale plot height so bars don't collapse with many samples
            bar_px  <- max(18L, round(500L / max(nrow(df), 1L) * nrow(df)))
            plot_h  <- paste0(max(300L, min(bar_px * nrow(df), 1400L)), "px")

            p <- ggplot(df, aes(x = library_size / 1e6, y = sample,
                                fill = condition, text = paste0(
                                    "Sample: ", sample,
                                    "\nLibrary size: ", round(library_size / 1e6, 2), "M",
                                    "\nCondition: ", condition
                                ))) +
                geom_col() +
                scale_fill_manual(values = pal) +
                labs(x = "Library size (millions)", y = NULL, fill = "Condition") +
                theme_bw(base_size = 11)
            ggplotly(p, tooltip = "text", height = max(300L, min(bar_px * nrow(df), 1400L)))
        })

        # Count Distribution: show the pre-rendered log2 CPM boxplot from the pipeline
        # (all samples, all conditions, generated by qc_analysis.R).
        # The VST matrix has ~78k genes including unexpressed ones (VST floor ≈ -0.38),
        # so re-computing VST quantiles in the dashboard gives misleading flat values.
        # The pipeline PNG uses log2(CPM+1) on detected genes — the correct metric.
        output$plot_count_dist <- renderImage({
            png_path <- normalizePath(
                file.path(data()$qc_basedir, "qc_count_distribution.png"),
                mustWork = FALSE
            )
            if (!file.exists(png_path)) {
                # fallback: look relative to RESULTS_DIR
                png_path <- normalizePath(
                    file.path(RESULTS_DIR, "qc", "qc_count_distribution.png"),
                    mustWork = FALSE
                )
            }
            if (file.exists(png_path)) {
                list(src = png_path, contentType = "image/png",
                     style = "max-width:100%; height:auto;",
                     alt  = "Count distribution (log2 CPM)")
            } else {
                list(src = "", alt = "Count distribution not available", style = "display:none")
            }
        }, deleteFile = FALSE)

        # PCA: reactive state for selected samples.
        # Initialized to all samples on first data load; updated only when Apply is clicked.
        pca_selected <- reactiveVal(NULL)

        observe({
            vst <- data()$vst
            ss  <- data()$samplesheet
            if (nrow(vst) == 0 || nrow(ss) == 0) return()
            all_samps <- intersect(ss$sample, colnames(vst))
            updateSelectizeInput(session, "pca_samples",
                                 choices = all_samps, selected = all_samps,
                                 server = TRUE)
            if (is.null(pca_selected())) pca_selected(all_samps)
            # Build dynamic color-by choices from available samplesheet columns
            col_choices <- c("Condition" = "condition")
            if ("norm_group" %in% colnames(ss)) col_choices <- c(col_choices, "Norm group" = "norm_group")
            if ("batch"      %in% colnames(ss)) col_choices <- c(col_choices, "Batch"      = "batch")
            updateSelectInput(session, "pca_color_by",
                              choices = col_choices, selected = "condition")
        })

        observeEvent(input$pca_apply, {
            req(input$pca_samples)
            pca_selected(input$pca_samples)
        })

        observeEvent(input$pca_reset, {
            vst <- data()$vst
            ss  <- data()$samplesheet
            all_samps <- intersect(ss$sample, colnames(vst))
            updateSelectizeInput(session, "pca_samples", selected = all_samps)
            pca_selected(all_samps)
        })

        output$plot_pca <- renderPlotly({
            vst <- data()$vst
            ss  <- data()$samplesheet
            sel <- pca_selected()
            if (is.null(sel) || nrow(vst) == 0 || nrow(ss) == 0)
                return(plot_ly() %>% layout(title = "No data"))

            sample_cols <- intersect(sel, colnames(vst))
            if (length(sample_cols) < 2)
                return(plot_ly() %>% layout(title = "Select at least 2 samples"))

            mat <- t(as.matrix(vst[, sample_cols, drop = FALSE]))
            pca <- prcomp(mat, scale. = FALSE)
            pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

            color_var <- input$pca_color_by %||% "condition"
            if (!color_var %in% colnames(ss)) color_var <- "condition"

            pca_df <- data.frame(
                sample = rownames(pca$x),
                PC1    = pca$x[, 1],
                PC2    = pca$x[, 2]
            )
            pca_df[[color_var]] <- ss[[color_var]][match(pca_df$sample, ss$sample)]

            groups <- unique(pca_df[[color_var]])
            pal    <- setNames(
                hcl.colors(length(groups), palette = "Dark 3", fixup = TRUE),
                groups
            )

            hover_text <- paste0(
                "Sample: ",    pca_df$sample,
                "\nCondition: ", ss$condition[match(pca_df$sample, ss$sample)],
                if ("norm_group" %in% colnames(ss))
                    paste0("\nNorm group: ", ss$norm_group[match(pca_df$sample, ss$sample)]) else "",
                if ("batch" %in% colnames(ss))
                    paste0("\nBatch: ", ss$batch[match(pca_df$sample, ss$sample)]) else ""
            )

            n_shown <- length(sample_cols)
            n_total <- length(intersect(ss$sample, colnames(vst)))
            subtitle <- if (n_shown < n_total)
                paste0(n_shown, " of ", n_total, " samples selected")
            else
                paste0("All ", n_total, " samples")

            plot_ly(
                pca_df,
                x         = ~PC1,
                y         = ~PC2,
                color     = pca_df[[color_var]],
                colors    = pal,
                text      = hover_text,
                hoverinfo = "text",
                type      = "scatter",
                mode      = "markers",
                marker    = list(size = 10)
            ) %>%
                layout(
                    xaxis  = list(title = paste0("PC1 (", pct[1], "%)")),
                    yaxis  = list(title = paste0("PC2 (", pct[2], "%)")),
                    title  = list(
                        text = paste0("PCA \u2014 VST (blind)<br><sup>", subtitle, "</sup>")
                    ),
                    legend = list(title = list(text = color_var))
                )
        })

        output$corr_heatmap <- renderImage({
            img_path <- normalizePath(data()$qc_heatmap_path, mustWork = FALSE)
            if (!file.exists(img_path)) {
                list(src = "", alt = "Correlation heatmap not available",
                     style = "display:none")
            } else {
                list(src = img_path, contentType = "image/png",
                     alt = "Sample correlation heatmap",
                     style = "max-width:100%; height:auto;")
            }
        }, deleteFile = FALSE)
    })
}
