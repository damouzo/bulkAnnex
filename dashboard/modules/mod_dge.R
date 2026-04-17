# modules/mod_dge.R — Differential Expression tab

mod_dge_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12, h3("Differential Expression"))
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
                    numericInput(ns("n_label"),  "Label top N genes", value = 20, min = 0, max = 100),
                    sliderInput(ns("label_size"), "Label font size",
                                min = 6, max = 16, value = 9, step = 1)
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Volcano",
                        plotlyOutput(ns("volcano_plot"), height = "500px")
                    ),
                    nav_panel("MA",
                        plotlyOutput(ns("ma_plot"), height = "500px")
                    ),
                    nav_panel("Results Table",
                        DTOutput(ns("results_table"))
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

        current_dge <- reactive({
            req(input$contrast)
            df <- data()$dge[[input$contrast]]
            if (is.null(df)) return(data.frame())
            df %>%
                filter(!is.na(padj)) %>%
                mutate(
                    direction = case_when(
                        padj < input$padj_cut & log2FoldChange >= input$lfc_cut  ~ "Up",
                        padj < input$padj_cut & log2FoldChange <= -input$lfc_cut ~ "Down",
                        TRUE ~ "NS"
                    )
                )
        })

        output$volcano_plot <- renderPlotly({
            df <- current_dge()
            if (nrow(df) == 0) return(plot_ly() %>% layout(title = "No data"))

            # Top N labels
            n_lab  <- min(input$n_label, nrow(df))
            top_ids <- df %>% filter(direction != "NS") %>%
                arrange(padj) %>% head(n_lab) %>% pull(gene_id)
            df$label <- ifelse(df$gene_id %in% top_ids, df$gene_name, "")

            dir_col <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "#bbbbbb")

            plot_ly(df,
                    x    = ~log2FoldChange,
                    y    = ~(-log10(padj)),
                    text = ~paste0("Gene: ", gene_name,
                                   "<br>log2FC: ", round(log2FoldChange, 3),
                                   "<br>padj: ",   signif(padj, 3)),
                    color       = ~direction,
                    colors      = dir_col,
                    hoverinfo   = "text",
                    type        = "scatter",
                    mode        = "markers",
                    marker      = list(size = 4, opacity = 0.7)) %>%
                add_text(text = ~label, textposition = "top center",
                         textfont = list(size = input$label_size), showlegend = FALSE) %>%
                layout(
                    xaxis = list(title = "Shrunken log2 fold change"),
                    yaxis = list(title = "-log10(padj)"),
                    title = paste0("Volcano: ", input$contrast),
                    shapes = list(
                        list(type = "line", x0 = input$lfc_cut, x1 = input$lfc_cut,
                             y0 = 0, y1 = 1, yref = "paper",
                             line = list(color = "grey", dash = "dash")),
                        list(type = "line", x0 = -input$lfc_cut, x1 = -input$lfc_cut,
                             y0 = 0, y1 = 1, yref = "paper",
                             line = list(color = "grey", dash = "dash")),
                        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
                             y0 = -log10(input$padj_cut), y1 = -log10(input$padj_cut),
                             line = list(color = "grey", dash = "dash"))
                    )
                )
        })

        output$ma_plot <- renderPlotly({
            df <- current_dge()
            if (nrow(df) == 0) return(plot_ly() %>% layout(title = "No data"))

            n_lab   <- min(input$n_label, nrow(df))
            top_ids <- df %>% filter(direction != "NS") %>%
                arrange(padj) %>% head(n_lab) %>% pull(gene_id)
            df$label <- ifelse(df$gene_id %in% top_ids, df$gene_name, "")

            plot_ly(df,
                    x    = ~log10(baseMean + 1),
                    y    = ~log2FoldChange,
                    text = ~paste0("Gene: ", gene_name,
                                   "<br>baseMean: ", round(baseMean, 1),
                                   "<br>log2FC: ",   round(log2FoldChange, 3),
                                   "<br>padj: ",     signif(padj, 3)),
                    color     = ~direction,
                    colors    = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "#bbbbbb"),
                    hoverinfo = "text",
                    type      = "scatter",
                    mode      = "markers",
                    marker    = list(size = 4, opacity = 0.7)) %>%
                add_text(text = ~label, textposition = "top center",
                         textfont = list(size = input$label_size), showlegend = FALSE) %>%
                layout(
                    xaxis = list(title = "log10(mean normalised count + 1)"),
                    yaxis = list(title = "Shrunken log2 fold change"),
                    title = paste0("MA: ", input$contrast)
                )
        })

        output$results_table <- renderDT({
            df <- current_dge()
            if (nrow(df) == 0) return(datatable(data.frame()))
            df_show <- df %>%
                mutate(
                    log2FoldChange = round(log2FoldChange, 3),
                    baseMean       = round(baseMean, 1),
                    pvalue         = signif(pvalue, 3),
                    padj           = signif(padj, 3)
                ) %>%
                select(gene_id, gene_name, baseMean, log2FoldChange, pvalue, padj, direction)

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
