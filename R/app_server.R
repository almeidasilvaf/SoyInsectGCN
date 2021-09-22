#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @importFrom shiny updateSelectizeInput reactive icon renderPlot
#' @importFrom shinydashboard valueBox renderValueBox
#' @importFrom ggplot2 aes_ ggplot guides labs
#' @importFrom ggnetwork theme_blank geom_edges geom_nodes
#' geom_nodelabel_repel unit
#' @importFrom DT datatable renderDataTable
#' @noRd
app_server <- function(input, output, session) {
    shiny::updateSelectizeInput(session, 'gene_id', 
                                choices = genes_modules$Genes, server = TRUE)
    
    #----Start value boxes----
    # Module
    mod <- shiny::reactive({
        m <- genes_modules[genes_modules$Genes == input$gene_id, 2]
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
        d <- scaled_degree$Scaled[scaled_degree$Gene == input$gene_id]
        d
    })
    output$degree <- renderValueBox({
        valueBox(paste0(deg()),
                 "Scaled degree", color="red", icon=shiny::icon("bar-chart-o"))
    })
    # Hub status
    hstatus <- reactive({
        status <- "No"
        if(input$gene_id %in% hubs$Gene) { status <- "Yes" }
        status
    })
    output$hub_status <- renderValueBox({
        valueBox(paste0(hstatus()),
                 "Hub status", color="green", icon=shiny::icon("star"))
    })
    ## End of value boxes
    
    #----Add DT of module enrichment
    enr_df <- reactive({
        m <- genes_modules[genes_modules$Genes == input$gene_id, 2]
        enrich <- mod_enrich[mod_enrich$Module == m, c(1,2,4)]
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
    
    #----Add network visualization----
    
    output$netviz_static <- shiny::renderPlot({
        n <- plotdata[plotdata$Modules == mod(), ]
        n$isInput <- ifelse(n$name == input$gene_id, TRUE, FALSE)
        p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, 
                                             xend = ~xend, yend = ~yend)) + 
            ggnetwork::geom_edges(color = "grey75", alpha = 0.5, 
                                  show.legend = FALSE) +
            ggnetwork::geom_nodes(aes_(size = ~Scaled10, alpha = ~Scaled4), 
                                  color = mod()) +
            ggplot2::labs(size="Scaled degree * 10") +
            ggplot2::guides(alpha = "none") +
            ggnetwork::geom_nodelabel_repel(
                ggplot2::aes_(label = ~name),
                color="black", box.padding = ggnetwork::unit(1, "lines"), 
                data = function(x) {
                    x[x$isInput, ]
                }, show.legend = FALSE) +
            ggnetwork::theme_blank()
        p
    })
}
