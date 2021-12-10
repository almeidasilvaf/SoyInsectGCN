#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @importFrom shiny updateSelectizeInput reactive icon renderPlot observe
#' @importFrom shinydashboard valueBox renderValueBox
#' @importFrom ggplot2 aes_ ggplot guides labs
#' @importFrom ggnetwork theme_blank geom_edges geom_nodes
#' geom_nodelabel_repel unit
#' @importFrom DT datatable renderDataTable
#' @noRd
app_server <- function(input, output, session) {
    
    #----Which taxon: 'insect' or 'nematode'?-----------------------------------
    taxon <- shiny::reactive({
        t <- input$taxon_id
        t
    })
    
    shiny::observeEvent(taxon(), {
        gene_opts <- as.character(
            genes_modules[genes_modules$taxon == taxon(), 1]
        )
        shiny::updateSelectizeInput(inputId = "gene_id", choices = gene_opts,
                                    server = TRUE) 
    })
    
    #----Start value boxes----
    # Module
    mod <- shiny::reactive({
        m <- genes_modules[genes_modules$taxon == taxon(), ]
        m <- m[m$Genes == input$gene_id, 2]
        m
    })
    output$module <- shinydashboard::renderValueBox({
        shinydashboard::valueBox(
            paste0(mod()),
            "Module", color="yellow", icon=icon("project-diagram")
            )
    })
    # Scaled degree
    deg <- shiny::reactive({
        d <- scaled_degree[scaled_degree$taxon == taxon(), ]
        d <- d$Scaled[d$Gene == input$gene_id]
        d
    })
    output$degree <- renderValueBox({
        valueBox(paste0(deg()),
                 "Scaled degree", color="red", 
                 icon=shiny::icon("chart-bar"))
    })
    # Hub status
    hub_vector <- reactive({
        h <- as.character(hubs[[taxon()]])
        h
    })
    hstatus <- reactive({
        status <- "No"
        if(input$gene_id %in% hub_vector()) { status <- "Yes" }
        status
    })
    output$hub_status <- renderValueBox({
        valueBox(paste0(hstatus()),
                 "Hub status", color="green", icon=shiny::icon("star"))
    })
    ## End of value boxes
    
    #----Add DT of module enrichment
    enr_df <- reactive({
        enrich <- mod_enrich[mod_enrich$taxon == taxon(), ]
        enrich <- enrich[enrich$Module == mod(), c(1,2,4)]
        colnames(enrich) <- c("Term", "P_adj", "Category")
        rownames(enrich) <- seq_len(nrow(enrich))
        enrich
    })
    output$mod_enrich <- DT::renderDataTable(
        #DT::datatable(enr_df()),
        enr_df(),
        selection = 'single',
        style = "bootstrap",
        rownames = FALSE,
        filter = 'top',
        options = list(lengthMenu = c(5, 10, 25, 50, 100), pageLength = 5)
        )
}
