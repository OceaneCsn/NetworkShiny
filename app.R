library(shiny)
library(shinythemes)
library(DT)
library(visNetwork)
library(stringr)
library(ggplot2)
setwd("./")
source("Visu.R")


# todo : ggplotly? Mettre les statistiques du graphe? ranking genes? choix du GRN ds menu déroulant

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Network visualisation"),
    hr(),
    fixedRow(
        column(
            width = 6,
            visNetworkOutput("network", height = "1000px"),
        ),
        column(
          width = 6,
          tabsetPanel(type = "tabs",
                      tabPanel("Ontologies", DT::dataTableOutput("Ontologies")),
                      tabPanel("Heatmap", plotOutput("heatmap"), br(), plotOutput("expression_plot"))
                      )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    load("./DataNetworkGenieCO2LowNitrate.RData")
    load("./normalized.count_At.RData")
    data$edges$value <- data$edges$weight
    data$nodes$group <- ifelse(data$nodes$group == 1, "Transcription Factor", "Target Gene")
    
    output$network <- renderVisNetwork({
        graph <- visNetwork(nodes = data$nodes, edges = data$edges) %>% 
          visNodes(borderWidth=0.5) %>% 
          visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
          visPhysics(solver = "forceAtlas2Based", timestep = 1, minVelocity=10, 
                     maxVelocity = 10, stabilization = F)%>%
          visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = TRUE)%>% 
          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
                      visGroups(groupname = "Transcription Factor", size = 28,
                                color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
                      visGroups(groupname = "Target Gene", color = "#77EEAA")
        graph
    })
    
    output$SelectedGene <- renderPrint({
        print(input$click)
    })
    
    #output$SelectedGene <- renderPrint({print(input$mynetwork_selectedBy)})
    

    output$Ontologies <- DT::renderDataTable({
        if(is.null(input$click)){
            data$nodes[c("description", "Ontology")]
        }
        else{
          neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                      data$edges[grepl(input$click[1], data$edges$to),]$from))
            data$nodes[c(input$click[1], neighboors),c("description", "Ontology")]
        }})
    
    output$expression_plot <- renderPlot({
        getExpression(input$click)
    })
    
    output$heatmap <- renderPlot({
      if(is.null(input$click)){
        heatmapPerso(normalized.count, conds = "all", genes = c("AT1G01150", "AT1G01590"))
      }
      else{
        neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                    data$edges[grepl(input$click[1], data$edges$to),]$from))
        heatmapPerso(normalized.count, conds = "all", genes = c(input$click[1], neighboors))
      }
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

