# modules/mod_about.R — About tab

mod_about_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(12, h3("About bulkAnnex"))
        ),
        fluidRow(
            column(8,
                card(
                    card_header("What is bulkAnnex?"),
                p("bulkAnnex is a Nextflow DSL2 pipeline for downstream bulk RNA-seq analysis.
                   It picks up where nf-core/rnaseq leaves off — starting from the Salmon
                   gene count matrix and running QC, normalisation, differential expression
                   (DESeq2), gene set enrichment (GO, KEGG, Reactome),
                   and this interactive dashboard. Reproducible, containerised, and runs
                   identically on a laptop or a thousand-core HPC cluster.")
                )
            ),
            column(4,
                card(
                    card_header("Links"),
                    tags$ul(
                        tags$li(tags$a(href = "https://github.com/damouzo/bulkAnnex",
                                       target = "_blank",
                                       icon("github"), " GitHub repository")),
                        tags$li(tags$a(href = "https://damouzo.github.io/",
                                       target = "_blank",
                                       icon("user"), " Portfolio"))
                    )
                ),
                card(
                    card_header("Software versions"),
                    uiOutput(ns("versions_table"))
                )
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
