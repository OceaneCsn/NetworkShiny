library(shiny)
library(shinythemes)
library(DT)
library(visNetwork)
library(stringr)
library(ggplot2)
#source("./DEFunctions.R")
setwd("./")
source("Visu.R")

getExpression <- function(gene, conds = "all", specie = "At"){
  # Plots the expression levels of a given gene, using the normized.count data provoded.
  # conditions are all the columns of the data by default, or can be specified
  # biological replicated should be identified by _n
  if(specie == "At") load("normalized.count_At.RData")
  if(specie == "Sl") load("normalized.count_Sl.RData")
  if (length(conds) ==1){
    conds = colnames(normalized.count)
  }else{conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))}
  df <- normalized.count[gene, conds]
  library(reshape2)
  d<- melt(df, silent=T)
  d$group = str_split_fixed(rownames(d), "_", 2)[,1]

  p <- ggplot(data = d, aes(x=group, y=value, fill=group)) + geom_dotplot(binaxis = "y", stackdir = "center") +
    scale_fill_discrete(name = "Conditions")+
    theme(strip.text.x = element_text(size = 26), plot.title = element_text(size=22, face="bold"),
          legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=20),
          axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 26, angle = 320, hjust = 0, colour = "grey50"),
          axis.title=element_text(size=17)) + ylab("Normalized counts") +
    ggtitle(paste("Normalized expression for ", gene)) + xlab("- C : elevated CO2 - c : ambiant CO2 - N : 10mM nitrate - n : 0.5mM nitrate - F : iron - f : iron starvation")
  return(p)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Network visualisation"),
    hr(),
    fluidRow(
        column(
            width = 6, height = 6,
            visNetworkOutput("network"),
            
        ),
        column(
          width = 6,
          div(h3("Gene expression accros all conditions"), align = "center"),
          plotOutput("heatmap")
          
        )
    ),
    
    hr(),
    fluidRow(
        column(
            width = 6,
            div(h3("You selected"), align = "center"),
            DT::dataTableOutput("Ontologies")
            
        ),
        column(
          width = 6,
          div(h3("Gene expression accros all conditions"), align = "center"),
          plotOutput("expression_plot")
          
        )
    ),
    
    hr(),
    fluidRow(
        
    )
    
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    load("./DataNetworkGenieCO2LowNitrate.RData")
    load("./normalized.count_At.RData")
    data$edges$value <- data$edges$weight
    
    
    output$network <- renderVisNetwork({
        graph <- visNetwork(nodes = data$nodes, edges = data$edges, height = "30%", width = "100%") %>%
            visEdges(smooth = FALSE, arrows = 'to') %>% visPhysics(solver = "forceAtlas2Based", timestep = 1, minVelocity=10, maxVelocity = 10, stabilization = F)%>%
            visOptions(selectedBy = "group", 
                       highlightNearest = TRUE, 
                       nodesIdSelection  = TRUE, collapse = TRUE)%>% visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}"
                       )
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
      neighboors <- unique( union(data$edges[grepl(input$click[1], data$edges$from),]$to,
                                  data$edges[grepl(input$click[1], data$edges$to),]$from))
      heatmapPerso(normalized.count, conds = "all", genes = c(input$click[1], neighboors))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

