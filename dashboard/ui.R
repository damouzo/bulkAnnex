# ui.R — bulkAnnex dashboard UI

suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(DT)
    library(plotly)
})

# Source modules
lapply(list.files("modules", pattern = "\\.R$", full.names = TRUE), source)

ui <- page_navbar(
    title = tags$span(
        tags$img(src = "bulkannex_logo.png", height = "30px", style = "margin-right: 8px;"),
        "bulkAnnex"
    ),
    theme = bs_theme(
        bootswatch  = "flatly",
        primary     = "#2c3e50",
        secondary   = "#18bc9c",
        base_font   = font_google("Inter"),
        code_font   = font_google("Source Code Pro")
    ),
    header = tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "bulkannex.css")
    ),
    window_title = "bulkAnnex",

    # ---- Tab 1: Overview ----------------------------------------------------
    nav_panel(
        title = "Overview",
        icon  = icon("home"),
        mod_overview_ui("overview")
    ),

    # ---- Tab 2: QC ----------------------------------------------------------
    nav_panel(
        title = "QC",
        icon  = icon("chart-bar"),
        mod_qc_ui("qc")
    ),

    # ---- Tab 3: Differential Expression ------------------------------------
    nav_panel(
        title = "Differential Expression",
        icon  = icon("dna"),
        mod_dge_ui("dge")
    ),

    # ---- Tab 4: GSEA --------------------------------------------------------
    nav_panel(
        title = "GSEA",
        icon  = icon("project-diagram"),
        mod_gsea_ui("gsea")
    ),

    # ---- Tab 5: Gene Explorer -----------------------------------------------
    nav_panel(
        title = "Gene Explorer",
        icon  = icon("search"),
        mod_gene_explorer_ui("gene_explorer")
    ),

    # ---- Tab 6: About -------------------------------------------------------
    nav_panel(
        title = "About",
        icon  = icon("info-circle"),
        mod_about_ui("about")
    ),

    # ---- Footer spacer -------------------------------------------------------
    nav_spacer(),
    nav_item(
        tags$small(
            class = "text-muted",
            paste0("bulkAnnex v1.0 | ", format(Sys.Date(), "%Y"))
        )
    )
)
