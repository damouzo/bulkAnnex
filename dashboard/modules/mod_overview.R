# modules/mod_overview.R — Overview tab

mod_overview_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12,
                h3("Experiment Overview"),
                p("Summary of the bulkAnnex analysis run.")
            )
        ),
        fluidRow(
            column(3,
                value_box(
                    title    = "Samples",
                    value    = textOutput(ns("n_samples")),
                    showcase = icon("vials"),
                    theme    = "primary"
                )
            ),
            column(3,
                value_box(
                    title    = "Conditions",
                    value    = textOutput(ns("n_conditions")),
                    showcase = icon("layer-group"),
                    theme    = "secondary"
                )
            ),
            column(3,
                value_box(
                    title    = "Contrasts",
                    value    = textOutput(ns("n_contrasts")),
                    showcase = icon("not-equal"),
                    theme    = "info"
                )
            ),
            column(3,
                value_box(
                    title    = "DEGs (any contrast)",
                    value    = textOutput(ns("n_degs")),
                    showcase = icon("dna"),
                    theme    = "success"
                )
            )
        ),
        fluidRow(
            column(6,
                card(
                    card_header("Samplesheet"),
                    DTOutput(ns("samplesheet_table"))
                )
            ),
            column(6,
                card(
                    card_header("Contrasts summary"),
                    DTOutput(ns("contrasts_summary"))
                )
            )
        )
    )
}

mod_overview_server <- function(id, app_data) {
    moduleServer(id, function(input, output, session) {
        data <- reactive(app_data())

        output$n_samples    <- renderText(nrow(data()$samplesheet))
        output$n_conditions <- renderText(length(unique(data()$samplesheet$condition)))
        output$n_contrasts  <- renderText(length(data()$dge))
        output$n_degs <- renderText({
            if (length(data()$dge) == 0) return(0)
            all_sig <- sapply(data()$dge, function(df) {
                if (is.null(df) || !"padj" %in% colnames(df)) return(0L)
                sum(df$padj < 0.05, na.rm = TRUE)
            })
            format(max(all_sig), big.mark = ",")
        })

        output$samplesheet_table <- renderDT({
            datatable(data()$samplesheet,
                      options = list(pageLength = 25, scrollX = TRUE),
                      rownames = FALSE, class = "compact stripe")
        })

        output$contrasts_summary <- renderDT({
            if (length(data()$dge) == 0) return(datatable(data.frame()))
            ct_meta <- data()$contrasts_meta
            has_ng  <- !is.null(ct_meta) && nrow(ct_meta) > 0 &&
                       "norm_group" %in% colnames(ct_meta) &&
                       "contrast_id" %in% colnames(ct_meta)
            summary_df <- do.call(rbind, lapply(names(data()$dge), function(cid) {
                df <- data()$dge[[cid]]
                ng <- if (has_ng) {
                    ct_meta$norm_group[ct_meta$contrast_id == cid][1] %||% "\u2014"
                } else "\u2014"
                if (is.null(df) || !"padj" %in% colnames(df)) {
                    return(data.frame(contrast = cid, norm_group = ng, up = 0, down = 0, total = 0))
                }
                data.frame(
                    contrast   = cid,
                    norm_group = ng,
                    up         = sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE),
                    down       = sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE),
                    total      = sum(df$padj < 0.05, na.rm = TRUE)
                )
            }))
            datatable(summary_df,
                      colnames = c("Contrast", "Norm group", "Up", "Down", "Total DEGs"),
                      options  = list(pageLength = 10),
                      rownames = FALSE, class = "compact stripe")
        })
    })
}
