library(shiny)
library(shinythemes)
library(DT)
library(visNetwork)
library(stringr)
library(ggplot2)
library(shinydashboard)
library(igraph)
setwd("./")
source("Visu.R")
#options(bitmapType='cairo')

# todo : ggplotly Mettre les statistiques du graphe? ranking genes?

# Define UI for application that draws a histogram
ui <- dashboardPage(skin="black",
  
    dashboardHeader(title = "Network visualisation"),
    
    dashboardSidebar(
      sidebarMenu(
        menuItem("Network", tabName = "Network", icon = icon("project-diagram")),
        menuItem("Expression database", tabName = "Expression_data", icon = icon("seedling"))
      )
    ),
    
    dashboardBody(
      
      tabItems(
        # First tab content
        tabItem(tabName = "Network",
          selectInput("select", label = h3("Select network"), width = 2000,
                      choices = list("Common resopnse to CO2 (131 genes)" = 1, "Response to CO2 in nitrate starvation (1100 genes)" = 2), 
                      selected = 1),
          hr(),

          fixedRow(
              column(
                  width = 6,
                  textInput("geneVis", "Visualize a precise gene (AGI)", value = ""),
                  visNetworkOutput("network", height = "1000px"),
              ),
              column(
                width = 6,
                tabsetPanel(type = "tabs",
                            tabPanel("Ontologies", DT::dataTableOutput("Ontologies")),
                            tabPanel("Heatmap", plotOutput("heatmap"), br(), plotOutput("expression_plot")),
                            tabPanel("Gene ranking", DT::dataTableOutput("Ranking")),
                            tabPanel("Network statistics", plotOutput("NetStats"))
                            )
              )
          )
        ),
        tabItem(tabName = "Expression_data",
              textInput("gene", "Ask me a gene! (AGI)", value = "AT1G01020"),
              br(),
              plotOutput("expression_plot_specific")
              
        )
      )
    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    #load("./DataNetworkGenieCO2LowNitrate.RData")
    load("./normalized.count_At.RData")
    
    
    
    output$network <- renderVisNetwork({
        if(input$select ==1){
          load("./DataNetworkGenieCO2Clusters.RData")
        }
        else{
          load("./DataNetworkGenieCO2LowNitrate.RData")
        }
        data$edges$value <- data$edges$weight
        data$nodes$group <- ifelse(data$nodes$group == 1, "Transcription Factor", "Target Gene")
        
        
        
        visNetwork(nodes = data$nodes, edges = data$edges)%>% 
          visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
          visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10, 
                     maxVelocity = 10, stabilization = F)%>%
          visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = TRUE)%>% 
          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
                      visGroups(groupname = "Transcription Factor", size = 28,
                                color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
                      visGroups(groupname = "Target Gene", color = "#77EEAA") %>% 
          visNodes(borderWidth=0.5) 

    })
    
    observe({
      #input$click <- input$geneVis
      visNetworkProxy("network") %>%
        visFocus(id = input$geneVis, scale = 1) %>% visSelectNodes(id = input$geneVis)
    })
    
    output$SelectedGene <- renderPrint({
        print(input$click)
    })

    output$Ontologies <- DT::renderDataTable({
      if(input$select ==1){
        load("./DataNetworkGenieCO2Clusters.RData")
      }
      else{
        load("./DataNetworkGenieCO2LowNitrate.RData")
      }
        if(is.null(input$click)){
            data$nodes[c("description", "Ontology")]
        }
        else{
          neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                      data$edges[grepl(input$click[1], data$edges$to),]$from))
            data$nodes[c(input$click[1], neighboors),c("description", "Ontology")]
        }})
    
    output$Ranking <- DT::renderDataTable({
      if(input$select ==1){
        load("./DataNetworkGenieCO2Clusters.RData")
      }
      else{
        load("./DataNetworkGenieCO2LowNitrate.RData")
      }
      net <- igraph::graph_from_data_frame(d = data$edges)
      importance <- (degree(net)/max(degree(net))+betweenness(net)/max(betweenness(net)))*0.5
      data$nodes$Betweenness <- betweenness(net)[match(data$nodes$id, names(betweenness(net)))]
      data$nodes$Degree <- degree(net)[match(data$nodes$id, names(degree(net)))]
      data$nodes$Centrality <- importance[match(data$nodes$id, names(importance))]
      data$nodes$Rank <- order(-data$nodes$Centrality)
      
      data$nodes[order(-data$nodes$Centrality),c("description", "Ontology", "Centrality", "Betweenness", "Degree")]
      

    })
    
    output$expression_plot <- renderPlot({
        getExpression(input$click)
      
    })
    
    output$NetStats <- renderPlot({
      if(input$select ==1){
        load("./DataNetworkGenieCO2Clusters.RData")
      }
      else{
        load("./DataNetworkGenieCO2LowNitrate.RData")
      }
      net <- igraph::graph_from_data_frame(d = data$edges)
      netStats(net)
    })
    
    output$heatmap <- renderPlot({
      if(input$select ==1){
        load("./DataNetworkGenieCO2Clusters.RData")
      }
      else{
        load("./DataNetworkGenieCO2LowNitrate.RData")
      }
      if(is.null(input$click)){
        heatmapPerso(normalized.count, conds = "all", genes = c("AT1G01150", "AT1G01590"))
      }
      else{
        neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                    data$edges[grepl(input$click[1], data$edges$to),]$from))
        heatmapPerso(normalized.count, conds = "all", genes = c(input$click[1], neighboors))
      }
      
    })
    
    output$expression_plot <- renderPlot({
      getExpression(input$click)
    })
    
    output$expression_plot_specific <- renderPlot({
      getExpression(input$gene)
    })
    
    output$heatmap <- renderPlot({
      if(input$select ==1){
        load("./DataNetworkGenieCO2Clusters.RData")
      }
      else{
        load("./DataNetworkGenieCO2LowNitrate.RData")
      }
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

