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
                        plotlyOutput(ns("plot_library_sizes"), height = "400px")
                    ),
                    nav_panel("Count Distribution",
                        plotlyOutput(ns("plot_count_dist"), height = "400px")
                    ),
                    nav_panel("PCA",
                        plotlyOutput(ns("plot_pca"), height = "500px")
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
            if (nrow(df) == 0) return(datatable(data.frame()))
            df$library_size <- format(df$library_size, big.mark = ",", scientific = FALSE)
            df$n_detected   <- format(df$n_detected,   big.mark = ",", scientific = FALSE)
            datatable(df, options = list(pageLength = 20, scrollX = TRUE),
                      rownames = FALSE, class = "compact stripe")
        })

        output$plot_library_sizes <- renderPlotly({
            df <- qc_metrics()
            if (nrow(df) == 0) return(plot_ly() %>% layout(title = "No data"))
            df$sample <- factor(df$sample, levels = df$sample[order(df$library_size)])
            p <- ggplot(df, aes(x = library_size / 1e6, y = sample,
                                fill = condition, text = paste0(
                                    "Sample: ", sample,
                                    "\nLibrary size: ", round(library_size / 1e6, 2), "M",
                                    "\nCondition: ", condition
                                ))) +
                geom_col() +
                scale_fill_brewer(palette = "Set1") +
                labs(x = "Library size (millions)", y = NULL, fill = "Condition") +
                theme_bw(base_size = 11)
            ggplotly(p, tooltip = "text")
        })

        output$plot_count_dist <- renderPlotly({
            vst <- data()$vst
            ss  <- data()$samplesheet
            if (nrow(vst) == 0 || nrow(ss) == 0) return(plot_ly() %>% layout(title = "No data"))

            sample_cols <- ss$sample
            mat <- as.matrix(vst[, sample_cols, drop = FALSE])
            q_df <- data.frame(
                sample    = rep(sample_cols, each = 5),
                q         = rep(c(0, 0.25, 0.5, 0.75, 1), times = length(sample_cols)),
                value     = as.vector(apply(mat, 2, quantile, probs = c(0, .25, .5, .75, 1)))
            )
            q_df$condition <- ss$condition[match(q_df$sample, ss$sample)]
            p <- plot_ly(data = q_df %>% filter(q == 0.5),
                         x = ~sample, y = ~value, type = "bar",
                         color = ~condition,
                         text = ~paste0("Median VST: ", round(value, 2)),
                         hoverinfo = "text") %>%
                layout(title = "Median VST per sample", xaxis = list(tickangle = -45))
            p
        })

        output$plot_pca <- renderPlotly({
            vst <- data()$vst
            ss  <- data()$samplesheet
            if (nrow(vst) == 0 || nrow(ss) == 0) return(plot_ly() %>% layout(title = "No data"))

            sample_cols <- ss$sample
            mat <- t(as.matrix(vst[, sample_cols, drop = FALSE]))
            pca <- prcomp(mat, scale. = FALSE)
            pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
            pca_df <- data.frame(
                sample    = rownames(pca$x),
                PC1       = pca$x[, 1],
                PC2       = pca$x[, 2],
                condition = ss$condition[match(rownames(pca$x), ss$sample)]
            )
            has_batch <- "batch" %in% colnames(ss)
            if (has_batch) pca_df$batch <- ss$batch[match(pca_df$sample, ss$sample)]

            hover_text <- paste0("Sample: ", pca_df$sample,
                                 "\nCondition: ", pca_df$condition,
                                 if (has_batch) paste0("\nBatch: ", pca_df$batch) else "")

            plot_ly(pca_df, x = ~PC1, y = ~PC2,
                    color = ~condition,
                    text = hover_text,
                    hoverinfo = "text",
                    type = "scatter", mode = "markers+text",
                    textposition = "top center",
                    marker = list(size = 10)) %>%
                layout(
                    xaxis = list(title = paste0("PC1 (", pct[1], "%)")),
                    yaxis = list(title = paste0("PC2 (", pct[2], "%)")),
                    title = "PCA — VST (blind)"
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
