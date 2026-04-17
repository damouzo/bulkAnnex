# server.R — bulkAnnex dashboard server
# Libraries and data-loading helpers are defined in global.R (sourced by app.R).

# Source all module files
lapply(list.files("modules", pattern = "\\.R$", full.names = TRUE), source)

server <- function(input, output, session) {
    # Load all pipeline results once at startup; RESULTS_DIR comes from global.R
    app_data <- reactive({
        load_all_data(RESULTS_DIR)
    })

    # ---- Module servers -----------------------------------------------------
    mod_overview_server("overview",           app_data)
    mod_qc_server("qc",                       app_data)
    mod_dge_server("dge",                     app_data)
    mod_gsea_server("gsea",                   app_data)
    mod_gene_explorer_server("gene_explorer", app_data)
    mod_about_server("about",                 app_data)
}
