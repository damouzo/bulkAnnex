#!/usr/bin/env Rscript
# app.R — bulkAnnex Shiny dashboard entry point

source("global.R")
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)
