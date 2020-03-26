library(shiny)
library(shinythemes)
library(DT)
library(visNetwork)
library(stringr)
library(ggplot2)
library(shinydashboard)
library(igraph)
library(plotly)
if(version$os=="linux-gnu") options(bitmapType='cairo')


# personnal functions
source("./Functions/Visu.R")
source("./Functions/Network_functions.R")

# expression levels and ontologies
load("./Data/OntologyAllGenes.RData")
load("./Data/normalized.count_At.RData")
load("./Data/NitrateGenes.RData")
load("./Data/DEGsListsFiltered.RData")

listFiles <- list.files("./NetworkData/", full.names = F)
names(listFiles) = listFiles
files <- lapply(split(listFiles, names(listFiles)), unname)

load(paste0("./NetworkData/NitrateDEGenes_CO2-N.RData"))

load("./NetworkReference/Gaudinier_Nature_2018_TFs_Nitrates_Ref.RData")


ui <- dashboardPage(skin="black",
  
    dashboardHeader(title = "Network visualisation"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Network", tabName = "Network", icon = icon("project-diagram")),
        menuItem("Expression database", tabName = "Expression_data", icon = icon("seedling")),
        menuItem("Differential Expression database", tabName = "DEGs", icon = icon("seedling"))
      )
    ),
    
    dashboardBody(
      
      tabItems(
        tabItem(tabName = "Network",
                selectInput("select", label = h3("Select network data"), width = 700,
                            choices = files, selected="NitrateDEGenes_CO2-N.RData"),
                
          hr(),
          fixedRow(
              column(
                  width = 6,
                  textInput("geneVis", "Visualize a precise gene (AGI)", value = ""),
                  visNetworkOutput("network", height = "1000px"),
                  actionButton("colorFromNitrate", "Color according to nitrate score")
              ),
              column(
                width = 6,
                tabsetPanel(type = "tabs",
                            tabPanel("Ontologies", DT::dataTableOutput("Ontologies")),
                            tabPanel("Ontologies Nitrate", DT::dataTableOutput("nitrate")),
                            tabPanel("Heatmap", plotlyOutput("heatmap"), br(), plotOutput("expression_plot")),
                            tabPanel("Gene ranking", DT::dataTableOutput("Ranking")),
                            tabPanel("Network statistics", plotOutput("NetStats")),
                            tabPanel("Comparison to gaudinier",
                                     h2("Combinatorial interactions between transcription factors and promoters of genes associated with nitrogen metabolism, signalling and nitrogen-associated processes"),
                                     h3("Common links with state of the art network :"),
                                     DT::dataTableOutput("commonLinks"), h3('Common genes with state of the art network : '),
                                     DT::dataTableOutput("commonNodes"))
                            
                            #Pearson and Spearman Rank Correlation of transcription factor and target interactions across publically available nitrogen availability microarray experiments
                            )
              )
          )
        ),
        tabItem(tabName = "Expression_data",
              textInput("gene", "Ask me a gene! (AGI)", value = "AT1G01020"),
              br(),
              plotOutput("expression_plot_specific")
              
        ),
        
        tabItem(tabName = "DEGs",
                fixedRow(
                  column(
                    width = 10,
                    selectInput("comparison", label=h3("Select differentiel expression analysis"), choices = names(DEGs)),
                    DT::dataTableOutput("genesComp"))
                  ),
                  column(
                    width = 4,
                    textInput("geneToSearch", "Ask me a gene! (AGI)", value = "AT1G08090"),
                    verbatimTextOutput("comparisonsDEG")
                  )
                )
      
   
        )
      )
    )




# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
    #reactive({load(paste0("./NetworkData/",input$select))})
    
    output$network <- renderVisNetwork({
      load(paste0("./NetworkData/",input$select))
      
      
      
      data$nodes$label <- data$nodes$Ontology
      visNetwork(nodes = data$nodes, edges = data$edges)%>% 
        visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
        # visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10, 
        #            maxVelocity = 10, stabilization = F)%>%
        visPhysics(solver = "forceAtlas2Based", timestep = 0.6, minVelocity=12, 
                   maxVelocity = 10, stabilization = F)%>% 
        visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = F)%>% 
        visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
        visGroups(groupname = "Regulator", size = 28,
                  color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
        visGroups(groupname = "Target Gene", color = list("background" = "#77EEAA", hover = "grey")) %>% 
        visNodes(borderWidth=0.5, font = list("size"=35))  %>% visInteraction(hover = TRUE)

    })
    
    observe({
      visNetworkProxy("network") %>%
        visFocus(id = input$geneVis, scale = 1) %>% visSelectNodes(id = input$geneVis)
    })
    
    output$SelectedGene <- renderPrint({
        print(input$click)
    })

    output$Ontologies <- DT::renderDataTable({
      load(paste0("./NetworkData/",input$select))
      

        if(is.null(input$click)){
            data$nodes
        }
        else{
          neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                      data$edges[grepl(input$click[1], data$edges$to),]$from))
            data$nodes[c(input$click[1], neighboors),]
        }})
    
    output$nitrate <- DT::renderDataTable({
      load(paste0("./NetworkData/",input$select))
      
      
      if(is.null(input$click)){
        res <- NitrateGenes(data$nodes$id, nGenes, ontologies)
        res[,!grepl("_", colnames(res))]
      }
      else{
        neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                    data$edges[grepl(input$click[1], data$edges$to),]$from))
        res <- NitrateGenes(c(input$click[1], neighboors), nGenes, ontologies)
        res[,!grepl("_", colnames(res))]
      }})
    
    output$Ranking <- DT::renderDataTable({
      load(paste0("./NetworkData/",input$select))
      net <- igraph::graph_from_data_frame(d = data$edges)
      importance <- (degree(net)/max(degree(net))+betweenness(net)/max(betweenness(net)))*0.5
      data$nodes$Betweenness <- betweenness(net)[match(data$nodes$id, names(betweenness(net)))]
      data$nodes$Degree <- degree(net)[match(data$nodes$id, names(degree(net)))]
      data$nodes$Centrality <- importance[match(data$nodes$id, names(importance))]
      data$nodes$Rank <- order(-data$nodes$Centrality)
      data$nodes[order(-data$nodes$Centrality),c("description", "Ontology", "Centrality", "Betweenness", "Degree")]
    
    })
    
    output$expression_plot <- renderPlot({
        getExpression(input$click, normalized.count = normalized.count)
      
    })
    
    
    output$NetStats <- renderPlot({
      load(paste0("./NetworkData/",input$select))
      net <- igraph::graph_from_data_frame(d = data$edges)
      netStats(net)
    })
    
    
    observeEvent(input$colorFromNitrate,{
      load(paste0("./NetworkData/",input$select))
      data$nodes$label <- data$nodes$Ontology
      newNodes <- data$nodes
      nScores <- NitrateGenes(data$nodes$id, nGenes, ontologies)
      newNodes$group<- nScores[match(newNodes$id, nScores$Gene), "NitrateScore"]
      visNetworkProxy("network") %>%
        visUpdateNodes(nodes = newNodes) %>% visGroups(groupname = "6", size = 28,
                                                       color = "#660000") %>% 
        visGroups(groupname = "5", color = "#990000")%>% 
      visGroups(groupname = "4", color = "#cc0000")%>% 
      visGroups(groupname = "3", color = "#e06666")%>% 
      visGroups(groupname = "2", color = "#ea9999")%>% 
        visGroups(groupname = "1", color = "#f4cccc")%>% 
        visGroups(groupname = "0", color = "#ffffff")%>% 
        visGroups(groupname = "7", color = "#330000")
    })
    
    
    output$heatmap <- renderPlotly({
      load(paste0("./NetworkData/",input$select))
      if(is.null(input$click)){
        heatmapPerso(normalized.count, conds = "all", genes = c("AT1G01150", "AT1G01590"))
      }
      else{
        neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                    data$edges[grepl(input$click[1], data$edges$to),]$from))
        heatmapPerso(normalized.count, conds = "all", genes = c(input$click[1], neighboors))
      }
      
    })
    
    ################################ Functions DEG database
  
    
    output$expression_plot_specific <- renderPlot({
      getExpression(input$gene, normalized.count)
    })
    
    ################################ Functions DEG database
    
    output$genesComp <- DT::renderDataTable({
      
      NitrateGenes(DEGs[[input$comparison]], nGenes, ontologies)
    })
    
    output$comparisonsDEG <- renderPrint({
      res <- ""
      for(comp in names(DEGs)){
        if (input$geneToSearch %in% DEGs[[comp]])
          print(comp)
      }
      
    })
    
    output$commonLinks <- DT::renderDataTable({
      load(paste0("./NetworkData/",input$select))
      data$edges$pairs <- paste(data$edges$from, data$edges$to)
      
      
      gaudinier$edges$pairs <- paste(gaudinier$edges$from, gaudinier$edges$to)
      
      commonLinks <- intersect(gaudinier$edges$pairs, data$edges$pairs)
      data$edges[data$edges$pairs %in% commonLinks,c("from", "to")]
    })
    
    output$commonNodes <- DT::renderDataTable({
      load(paste0("./NetworkData/",input$select))
      data$nodes[data$nodes$id %in% gaudinier$nodes$id,c("Ontology", "description", "group", "ranking")]
    })
}


# Run the application 
shinyApp(ui = ui, server = server)