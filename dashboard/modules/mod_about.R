# modules/mod_about.R — About tab

mod_about_ui <- function(id) {
    ns <- NS(id)

    fluidRow(
        column(8,
            card(
                card_header("About bulkAnnex"),
                HTML(paste0(
                    "<div style='display:flex; align-items:center; gap:16px; margin-bottom:16px;'>",
                    "<img src='bulkAnnex_logo.png' height='56px'>",
                    "<div><h4 style='margin:0;'>bulkAnnex v1.0</h4>",
                    "<p style='margin:0; color:#666;'>Bulk RNA-seq downstream analysis pipeline &amp; dashboard</p></div>",
                    "</div>",
                    "<p>Interactive visualisation and exploration dashboard for bulk RNA-seq data ",
                    "processed through the <strong>bulkAnnex</strong> Nextflow DSL2 pipeline. ",
                    "Picks up where nf-core/rnaseq leaves off — starting from the Salmon gene ",
                    "count matrix and running normalisation, differential expression, and gene ",
                    "set enrichment. Reproducible, containerised, and runs identically on a ",
                    "laptop or an HPC cluster.</p>",
                    "<h5>Pipeline modules</h5>",
                    "<ul>",
                    "<li>Input validation (samplesheet, counts matrix, contrasts)</li>",
                    "<li>Quality control: library sizes, count distribution, PCA, correlation heatmap</li>",
                    "<li>Normalisation per norm_group (DESeq2 VST)</li>",
                    "<li>Differential expression (DESeq2 + lfcShrink): volcano, MA plot, heatmap</li>",
                    "<li>Gene set enrichment (fgsea / clusterProfiler): GO BP/MF/CC, KEGG, Reactome</li>",
                    "</ul>",
                    "<h5>Dashboard tabs</h5>",
                    "<ul>",
                    "<li><strong>Overview:</strong> Sample table, norm_group summary, contrast overview</li>",
                    "<li><strong>QC:</strong> Library sizes, count distribution, PCA, correlation heatmap</li>",
                    "<li><strong>DGE:</strong> Volcano plots with interactive labelling, results table</li>",
                    "<li><strong>GSEA:</strong> Dotplot, ridgeplot, running score, results table</li>",
                    "<li><strong>Gene Explorer:</strong> Per-gene expression across samples and conditions</li>",
                    "</ul>",
                    "<hr>",
                    "<p>",
                    "<a href='https://github.com/damouzo/bulkAnnex' target='_blank' style='margin-right:16px;'>",
                    "<i class='fab fa-github'></i> GitHub Repository</a>",
                    "</p>",
                    "<p class='text-muted' style='font-size:0.85rem;'>",
                    "<em>bulkAnnex v1.0 &nbsp;|&nbsp; Nextflow DSL2 &nbsp;+&nbsp; R/DESeq2/clusterProfiler &nbsp;+&nbsp; R/Shiny/bslib</em></p>"
                ))
            )
        ),
        column(4,
            card(
                card_header("Software versions"),
                uiOutput(ns("versions_table"))
            )
        )
    )
}

mod_about_server <- function(id, app_data) {
    moduleServer(id, function(input, output, session) {

        output$versions_table <- renderUI({
            # Load versions.yml files from the results directory
            ver_files <- list.files(RESULTS_DIR, pattern = "versions\\.yml$",
                                    full.names = TRUE, recursive = TRUE)

            # Parse YAML-ish lines: "    key: value"
            all_versions <- list()
            for (f in ver_files) {
                lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
                for (ln in lines) {
                    m <- regmatches(ln, regexec("^\\s+(\\S+):\\s+(.+)$", ln))[[1]]
                    if (length(m) == 3) all_versions[[m[2]]] <- m[3]
                }
            }

            # Known versions from the container (rocker/bioconductor:3.20)
            known <- list(
                "Nextflow"          = "24.10.x",
                "R"                 = "4.4.2",
                "Bioconductor"      = "3.20",
                "DESeq2"            = "1.46.x",
                "clusterProfiler"   = "4.14.x",
                "ReactomePA"        = "1.46.x",
                "pathview"          = "1.46.x",
                "enrichplot"        = "1.26.x",
                "ggridges"          = "0.5.x",
                "AnnotationDbi"     = "1.68.x",
                "shiny"             = "1.13.0",
                "bslib"             = "0.9.x",
                "DT"                = "0.33.x",
                "plotly"            = "4.10.x"
            )

            # Override known with parsed YAML versions
            for (k in names(all_versions)) {
                known[[k]] <- all_versions[[k]]
            }

            rows <- lapply(names(known), function(pkg) {
                tags$tr(
                    tags$td(tags$b(pkg), style = "padding: 3px 10px 3px 0;"),
                    tags$td(known[[pkg]], style = "color: #555; padding: 3px 0;")
                )
            })

            tags$table(style = "font-size: 0.88em; width: 100%;",
                       do.call(tagList, rows))
        })
    })
}
