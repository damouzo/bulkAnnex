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
                    numericInput(ns("n_label"),  "Label top N genes", value = 30, min = 0, max = 100),
                    sliderInput(ns("label_size"), "Label font size",
                                min = 6, max = 16, value = 9, step = 1),
                    actionButton(ns("apply"), "Apply settings",
                                 icon = icon("play"),
                                 class = "btn-primary w-100 mt-2"),
                    downloadButton(ns("download_volcano"), "Download volcano PNG",
                                   class = "btn-outline-secondary w-100 mt-2")
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Volcano",
                        imageOutput(ns("volcano_plot"), width = "100%", height = "700px")
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

        # Frozen params — only updates when Apply is clicked OR contrast changes.
        # This prevents every keystroke in numericInputs from triggering recomputation.
        params <- reactive({
            list(
                contrast   = input$contrast,
                padj_cut   = input$padj_cut,
                lfc_cut    = input$lfc_cut,
                n_label    = input$n_label,
                label_size = input$label_size
            )
        }) |> bindEvent(input$apply, input$contrast, ignoreNULL = TRUE, ignoreInit = FALSE)

        current_dge <- reactive({
            p <- params()
            req(p$contrast)
            df <- data()$dge[[p$contrast]]
            if (is.null(df)) return(data.frame())
            df %>%
                mutate(
                    direction = case_when(
                        !is.na(padj) & !is.na(log2FoldChange) &
                            padj < p$padj_cut & log2FoldChange >= p$lfc_cut  ~ "Up",
                        !is.na(padj) & !is.na(log2FoldChange) &
                            padj < p$padj_cut & log2FoldChange <= -p$lfc_cut ~ "Down",
                        TRUE ~ "NS"
                    )
                )
        })

        build_volcano_plot <- function(df_raw, p) {
            req(nrow(df_raw) > 0)

            df <- df_raw %>%
                filter(!is.na(padj), !is.na(log2FoldChange)) %>%
                mutate(
                    padj_plot = pmax(padj, .Machine$double.xmin),
                    neglog10_padj = -log10(padj_plot),
                    direction = case_when(
                        padj < p$padj_cut & log2FoldChange >= p$lfc_cut  ~ "Up",
                        padj < p$padj_cut & log2FoldChange <= -p$lfc_cut ~ "Down",
                        TRUE ~ "NS"
                    ),
                    label_base = ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, gene_id)
                )

            if (nrow(df) == 0) {
                return(
                    ggplot() +
                        annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                                 label = "No DGE rows available for volcano plot") +
                        theme_void()
                )
            }

            n_lab  <- min(p$n_label, nrow(df))
            sig_df <- df %>% filter(direction != "NS")
            top_ids <- if (nrow(sig_df) > 0 && n_lab > 0) {
                max_lfc <- max(abs(sig_df$log2FoldChange), na.rm = TRUE)
                max_nlp <- max(sig_df$neglog10_padj, na.rm = TRUE)
                sig_df %>%
                    mutate(score = sqrt(
                        (abs(log2FoldChange) / pmax(max_lfc, 1e-9))^2 +
                        (neglog10_padj / pmax(max_nlp, 1e-9))^2
                    )) %>%
                    arrange(desc(score)) %>%
                    head(n_lab) %>%
                    pull(gene_id)
            } else character(0)
            df$label <- ifelse(df$gene_id %in% top_ids, df$label_base, NA_character_)

            dir_col <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "#bbbbbb")

            ggplot(df, aes(x = log2FoldChange, y = neglog10_padj, colour = direction)) +
                geom_point(size = 0.8, alpha = 0.7) +
                ggrepel::geom_text_repel(aes(label = label),
                                         size = p$label_size / 3,
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.2,
                                         seed = 42,
                                         na.rm = TRUE) +
                geom_vline(xintercept = c(-p$lfc_cut, p$lfc_cut), linetype = "dashed", colour = "grey40") +
                geom_hline(yintercept = -log10(pmax(p$padj_cut, .Machine$double.xmin)), linetype = "dashed", colour = "grey40") +
                scale_colour_manual(values = dir_col) +
                labs(title = paste0("Volcano: ", p$contrast),
                     x = "Shrunken log2 fold change",
                     y = "-log10(padj)",
                     colour = NULL) +
                theme_bw(base_size = 12) +
                theme(plot.margin = ggplot2::margin(t = 20, r = 10, b = 5, l = 5, unit = "pt")) +
                coord_cartesian(clip = "off")
        }

        # Keep the fast workflow: changing contrast shows pre-rendered pipeline PNG.
        # Recompute only when Apply is clicked for the current contrast.
        applied_volcano_params <- reactiveVal(NULL)

        observeEvent(input$contrast, {
            applied_volcano_params(NULL)
        }, ignoreInit = TRUE)

        observeEvent(input$apply, {
            applied_volcano_params(list(
                contrast = input$contrast,
                padj_cut = input$padj_cut,
                lfc_cut = input$lfc_cut,
                n_label = input$n_label,
                label_size = input$label_size
            ))
        })

        output$volcano_plot <- renderImage({
            contrast <- input$contrast
            req(contrast)

            p_applied <- applied_volcano_params()
            apply_active <- !is.null(p_applied) && identical(p_applied$contrast, contrast)

            # 1) Default view: precomputed volcano from pipeline (instant, no recompute)
            if (!apply_active) {
                png_path <- get_dge_volcano_png(RESULTS_DIR, contrast)
                if (!is.null(png_path)) {
                    return(list(
                        src = png_path,
                        contentType = "image/png",
                        style = "width:100%; height:700px; object-fit:contain; background:white;"
                    ))
                }
            }

            # 2) Apply clicked (or missing precomputed PNG): build custom volcano
            df_raw <- data()$dge[[contrast]]
            if (is.null(df_raw) || nrow(df_raw) == 0) {
                tmpfile <- tempfile(fileext = ".png")
                p_empty <- ggplot() +
                    annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                             label = "No DGE results available for this contrast") +
                    theme_void()
                ggplot2::ggsave(tmpfile, p_empty, width = 8, height = 8, dpi = 150, bg = "white")
                return(list(
                    src = tmpfile,
                    contentType = "image/png",
                    style = "width:100%; height:700px; object-fit:contain; background:white;"
                ))
            }

            p_use <- if (apply_active) p_applied else list(
                contrast = contrast,
                padj_cut = input$padj_cut,
                lfc_cut = input$lfc_cut,
                n_label = input$n_label,
                label_size = input$label_size
            )

            tmpfile <- tempfile(fileext = ".png")
            p_vol <- build_volcano_plot(df_raw, p_use)
            ggplot2::ggsave(tmpfile, p_vol, width = 8, height = 8, dpi = 150, bg = "white")
            list(
                src = tmpfile,
                contentType = "image/png",
                style = "width:100%; height:700px; object-fit:contain; background:white;"
            )
        }, deleteFile = FALSE)

        output$download_volcano <- downloadHandler(
            filename = function() {
                contrast <- input$contrast %||% "contrast"
                paste0(contrast, "_volcano_dashboard.png")
            },
            content = function(file) {
                contrast <- input$contrast
                req(contrast)
                df_raw <- data()$dge[[contrast]]
                req(!is.null(df_raw), nrow(df_raw) > 0)

                p_applied <- applied_volcano_params()
                p_use <- if (!is.null(p_applied) && identical(p_applied$contrast, contrast)) {
                    p_applied
                } else {
                    list(
                        contrast = contrast,
                        padj_cut = input$padj_cut,
                        lfc_cut = input$lfc_cut,
                        n_label = input$n_label,
                        label_size = input$label_size
                    )
                }

                p_vol <- build_volcano_plot(df_raw, p_use)
                ggplot2::ggsave(file, p_vol, width = 8, height = 8, dpi = 300, bg = "white")
            }
        )

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
