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
load("./Data/PlnTFBDRegulatorsList.RData")

listFiles <- list.files("./NetworkData/", full.names = F)
names(listFiles) = listFiles
files <- lapply(split(listFiles, names(listFiles)), unname)

load(paste0("./NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData"))

load("./NetworkReference/Gaudinier_Nature_2018_TFs_Nitrates_Ref.RData")
load("./NetworkReference/Gaudinier_Correlations_Litterature.RData")


universe <- ontologies$entrezgene_id


colorEdge <- function(tested, validated){
  if(tested & validated) return("red")
  if(tested & !validated) return("black")
  if(!tested & !validated) return("white")
}

ui <- dashboardPage(skin="black",
                    
  
    dashboardHeader(title = "Network visualisation"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Network", tabName = "Network", icon = icon("project-diagram")),
        menuItem("Expression database", tabName = "Expression_data", icon = icon("seedling")),
        menuItem("Differential Expression database", tabName = "DEGs", icon = icon("table")),
        menuItem("Compare 2 networks", tabName = "Comparaison", icon = icon("greater-than-equal")),
        menuItem("Communities discovery", tabName = "Clustering", icon = icon("circle-notch"))
      )
    ),
    
    dashboardBody(
      tags$head(tags$style(".shiny-progress {top: 50% !important;left: 50% !important;margin-top: -100px !important;margin-left: -250px !important; color: blue;font-size: 20px;font-style: italic;}")),
      tabItems(
        tabItem(tabName = "Network",
                selectInput("select", label = h3("Select network data"), width = 700,
                            choices = files, selected="CO2DEGenes_faibleNitrate_CO2-N.RData"),
                
          hr(),
          fixedRow(
              column(
                  width = 5,
                  textInput("geneVis", "Visualize a precise gene (AGI)", value = ""),
                  visNetworkOutput("network", height = "800px"),
                  actionButton("colorFromNitrate", "Color according to nitrate score")
              ),
              column(
                width = 4,
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
                                     DT::dataTableOutput("commonNodes")),
                            tabPanel("DAP-Seq validated interactions", checkboxInput("edgesDap", "Color edges from DAP-Seq", FALSE), 
                                     h3("To check if the predicted target's promoters have the regulator's binding motif"),
                                     h2("Interactions validated by DAP-Seq :"),
                                     DT::dataTableOutput("edgesListDap", width = 700),
                                     h3("Taux de validation : "),
                                     textOutput("dapStats"),
                                     checkboxInput("nlpTargets", "Color NLP7 DAP-Seq targets", FALSE)
                                     )
                            
                            #Pearson and Spearman Rank Correlation of transcription factor and target interactions across publically available nitrogen availability microarray experiments
                            )
              )
          )
        ),
        tabItem(tabName = "Expression_data",
              textInput("gene", "Ask me a gene! (AGI)", value = "AT1G01020"),
              br(),
              plotOutput("expression_plot_specific", width="1400px") 
              
        ),
        
        tabItem(tabName = "DEGs",
                fixedRow(
                  column(
                    width = 6,
                    h3("See which genes aree differentially expressed in one transcriptome comparison :"),
                    selectInput("comparison", label=h3("Select differentiel expression analysis"), choices = names(DEGs)),
                    DT::dataTableOutput("genesComp", width = "400px"))
                  ),
                  column(
                    width = 4,
                    h3("For a gene, find in which transcriptome comparison it is differentially expressed :"),
                    textInput("geneToSearch", "Ask me a gene! (AGI)", value = "AT1G08090"),
                    verbatimTextOutput("comparisonsDEG")
                  )
                ),
        
        tabItem(tabName = "Comparaison",
                
                selectInput("select1", label = h3("Select first network data"), width = 600,
                            choices = files, selected="NitratesDEGenes_CO2-N.RData"),
                
                selectInput("select2", label = h3("Select second network data"), width = 600,
                            choices = files, selected="NitratesDEGenes_N.RData"),
                
                fixedRow(
                  column(
                    width = 5,
                    h3("Networks intersection :"),visNetworkOutput("netIntersection", height = "1000px")
                    )
                  ,
                  column(
                    width = 5,
                    h2("Common links :"),
                    DT::dataTableOutput("commonLinks12"), h2('Common genes : '),
                    DT::dataTableOutput("commonNodes12")
                  )
                )
        ),
        tabItem(tabName = "Clustering",
                
                column(width = 3, 
                       selectInput("selectClusterNetwork", label = h3("Select network data"), width = 600,
                                   choices = files, selected="CO2DEGenes_faibleNitrate_CO2-N.RData")),
                column(width = 3, selectInput("clustType", label = h3("Select type of clustering"), width = 600,
                                              choices = list("Topological clustering (Louvain)" = "louvain", "Expression clustering (COSEQ Poisson mixture models)"="coseq"),
                                              selected = "louvain")),
                
                column(width = 3,
                       selectInput("Module", label = h3("Select module"), width = 600,
                                   choices = NA)),
               
                fixedRow(
                  column(
                    width = 4,
                    h3("Network modules :"),visNetworkOutput("netClustering", height = "1000px")
                  )
                  ,
                  column(
                    width = 5,
                    tabsetPanel(type = "tabs",
                                tabPanel("Genes list", h2("Genes in that Community :"),
                                         DT::dataTableOutput("communityList")),
                                tabPanel("GO enrichment", h2("Ontologies enriched in this cluster : "), plotlyOutput("GOEnrich", width = 900, height = 900)),
                                tabPanel("GO enrichment comparison", h2("Comparison of the communitites : "), plotOutput("GOEnrichComp", height=900))
                    ),
                    
                  )
                )
          )
        )
      )
    )




# Define server logic required to draw a histogram
server <- function(input, output, session) {
  

  data <- reactive({
    load(paste0("./NetworkData/",input$select))
    data$nodes$label <- data$nodes$Ontology
    data$edges$color <-'#333366'
    data$edges$Regulator_Name <- ontologies[match(data$edges$from, ontologies$ensembl_gene_id),]$external_gene_name
    data$edges$Target_Name <- ontologies[match(data$edges$to, ontologies$ensembl_gene_id),]$external_gene_name
    data
  })
  
    output$network <- renderVisNetwork({
      
      visNetwork(nodes = data()$nodes, edges = data()$edges)%>% 
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
      visNetworkProxy("network") %>% visFocus(id = input$geneVis, scale = 1) %>% visSelectNodes(id = input$geneVis)
    })
    
    output$SelectedGene <- renderPrint({
        print(input$click)
    })

    output$Ontologies <- DT::renderDataTable({

        if(is.null(input$click)){
            data()$nodes
        }
        else{
          neighboors <- unique( union(data()$edges[grepl(input$click[1], data()$edges$from),]$to,
                                      data()$edges[grepl(input$click[1], data()$edges$to),]$from))
            data()$nodes[c(input$click[1], neighboors),]
        }})
    
    output$nitrate <- DT::renderDataTable({

      if(is.null(input$click)){
        res <- NitrateGenes(data()$nodes$id, nGenes, ontologies)
        res[,!grepl("_", colnames(res))]
      }
      else{
        neighboors <- unique( union(data()$edges[grepl(input$click[1], data()$edges$from),]$to,
                                    data()$edges[grepl(input$click[1], data()$edges$to),]$from))
        res <- NitrateGenes(c(input$click[1], neighboors), nGenes, ontologies)
        res[,!grepl("_", colnames(res))]
      }})
    
    output$Ranking <- DT::renderDataTable({

      net <- igraph::graph_from_data_frame(d = data()$edges)
      importance <- (degree(net)/max(degree(net))+betweenness(net)/max(betweenness(net)))*0.5
      
      data <- data()
      
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
      net <- igraph::graph_from_data_frame(d = data()$edges)
      netStats(net)
    })
    
    
    observeEvent(input$colorFromNitrate,{

      newNodes <- data()$nodes
      nScores <- NitrateGenes(data()$nodes$id, nGenes, ontologies)
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

      if(is.null(input$click)){
        heatmapPerso(normalized.count, conds = "all", genes = c("AT1G01150", "AT1G01590"))
      }
      else{
        neighboors <- unique( union(data()$edges[grepl(input$click[1], data()$edges$from),]$to,
                                    data()$edges[grepl(input$click[1], data()$edges$to),]$from))
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
    
    output$netIntersection <- renderVisNetwork({
      load(paste0("./NetworkData/",input$select1))
      network1 <- igraph::graph_from_data_frame(d = data$edges)
      
      load(paste0("./NetworkData/",input$select2))
      network2 <- igraph::graph_from_data_frame(d = data$edges)
      
      inter <- igraph::intersection(network1, network2, keep.all.vertices=F)
      dataInter <- networkData(inter, ontologies, TF)
      dataInter$nodes$label <- dataInter$nodes$Ontology
      
      nScores <- NitrateGenes(dataInter$nodes$id, nGenes, ontologies)
      dataInter$nodes$group<- nScores[match(dataInter$nodes$id, nScores$Gene), "NitrateScore"]
      
      
      visNetwork(nodes = dataInter$nodes, edges = dataInter$edges)%>% 
      visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
      visPhysics(solver = "forceAtlas2Based", timestep = 0.6, minVelocity=12, 
                 maxVelocity = 10, stabilization = F)%>% 
      visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = F)%>% 
      visEvents(click = "function(nodes){
                Shiny.onInputChange('click', nodes.nodes);
                ;}")%>%  visNodes(borderWidth=0.5, font = list("size"=35))  %>% visInteraction(hover = TRUE)%>% visGroups(groupname = "6", size = 28, color = "#660000") %>% 
      visGroups(groupname = "5", color = "#990000")%>% 
      visGroups(groupname = "4", color = "#cc0000")%>% 
      visGroups(groupname = "3", color = "#e06666")%>% 
      visGroups(groupname = "2", color = "#ea9999")%>% 
      visGroups(groupname = "1", color = "#f4cccc")%>% 
      visGroups(groupname = "0", color = "#ffffff")%>% 
      visGroups(groupname = "7", color = "#330000")
      
      
    })
    
    output$commonNodes12 <- DT::renderDataTable({
      
      load(paste0("./NetworkData/",input$select1))
      data1 <- data

      load(paste0("./NetworkData/",input$select2))

      data1$nodes[data1$nodes$id %in% data$nodes$id,c("Ontology", "description", "group", "ranking")]
    })
    
    output$commonLinks12 <- DT::renderDataTable({
      load(paste0("./NetworkData/",input$select1))
      network1 <- igraph::graph_from_data_frame(d = data$edges)
      
      load(paste0("./NetworkData/",input$select2))
      network2 <- igraph::graph_from_data_frame(d = data$edges)
      
      inter <- igraph::intersection(network1, network2, keep.all.vertices=F)
      dataInter <- networkData(inter, ontologies, TF)
      
      dataInter$edges[,c("from", "to")]
    })
    
    #################### Clustering network
    
    
    output$netClustering <- renderVisNetwork({
      plotNetwork(dataClust())
    })

    
    dataClust <- reactive({
      load(paste0("./NetworkData/",input$selectClusterNetwork))
      
      if(input$clustType == "louvain"){
        net <- igraph::graph_from_data_frame(d = data$edges, directed = F)
        communities <- cluster_louvain(net)
        membership <- membership(communities)
        data$nodes$group <- membership[match(data$nodes$id, names(membership))]
      }
      if(input$clustType == "coseq"){
        data$nodes$group <- data$nodes$coseqCluster
      }
      updateSelectInput(session, "Module", choices = data$nodes$group)
      data$nodes$label <- data$nodes$Ontology
      data$edges$color <-'#333366'
      data
    })
    
    observe({
      genes <- dataClust()$nodes[dataClust()$nodes$group ==input$Module,]$id
      visNetworkProxy("netClustering") %>% visSelectNodes( genes, highlightEdges = TRUE, clickEvent = T)%>% visOptions( highlightNearest = list("hideColor" = rgb(0.8,0.8,0.8,0.5)))
    })
    
    output$communityList <- DT::renderDataTable({
      dataClust()$nodes[dataClust()$nodes$group ==input$Module,]
    })
    
    output$GOEnrich <- renderPlotly({
      ids <- as.character(ontologies[match(dataClust()$nodes[dataClust()$nodes$group==input$Module,]$id, ontologies$ensembl_gene_id),]$entrezgene_id)
      print(ids)
      withProgress(message = "Ontologies Enrichment", {simOnt <- OntologyEnrich(ids, as.character(universe))})
      simOnt@result <- simOnt@result[order(-simOnt@result$p.adjust),]
      values <- str_split_fixed(simOnt@result$GeneRatio, "/", 2)
      simOnt@result$GeneRatio <- as.numeric(values[,1])/ as.numeric(values[,2])
      ggplotly(ggplot(data= simOnt@result, aes(x=Description, y=GeneRatio, fill=p.adjust)) + geom_bar(stat = "identity")+ coord_flip() +
                 ggtitle(paste("Enriched Ontologies for module", input$Module)) + ylab("") + xlab("Gene Ratio")+
                 theme(plot.title = element_text(size=15, face="bold"),
                       legend.title = element_text(size = 20, face="bold"), legend.text = element_text(size=15),
                       axis.text.y = element_text(size = 15, angle = 20), axis.text.x = element_text(size =8, angle = 0, hjust = 0, colour = "grey50"),
                       axis.title.y=element_blank()) )
    })
    
    output$GOEnrichComp <- renderPlot({
      idsList <- list()
      for(k in unique(dataClust()$nodes$group)){
        idsList[[as.character(k)]] <- na.omit(as.character(ontologies[match(dataClust()$nodes[dataClust()$nodes$group==k,]$id, ontologies$ensembl_gene_id),]$entrezgene_id))
      }
      withProgress(message = 'Ontologies enrichment comparison', {compareOnt(idsList=idsList, as.character(universe))})
    })
    
    #################################### DAPSeq functions
    
    observeEvent(input$edgesDap, {
      
      if(input$edgesDap){
        data <- data()
        data$edges$color <- mapply(colorEdge,  data$edges$testedByDapSeq, data$edges$DapSeqAproved)
      visNetworkProxy("network") %>% visUpdateEdges(data$edges)%>% 
        visGroups(groupname = "Regulator", size = 28,
                  color = list("background" = "grey", "border"="#CCCCCC"), shape = "square") %>% 
        visGroups(groupname = "Target Gene", color = "lightgrey") %>% 
        visNodes(borderWidth=0.5, font=list("size"=36)) 
       
      }
      
      else{
        visNetworkProxy("network")  %>%  visUpdateEdges(data()$edges) %>% visGroups(groupname = "Regulator", size = 28,
                                             color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
          visGroups(groupname = "Target Gene", color = "#77EEAA") 
      }
       
    })
    
    output$edgesListDap <- DT::renderDataTable({
      dataDap <-data()$edges[data()$edges$DapSeqAproved,]
      if(is.null(input$click)){
        dataDap[,c("from", "Regulator_Name", "to", "Target_Name")]
      }
      else{
       
        dataDap[grepl(input$click, dataDap$from),c("from", "Regulator_Name", "to", "Target_Name")]
      }
      
    })
    
    output$dapStats <- renderText({
      paste("Interactions comportant un TF etudie par DAP-Seq : ", round(sum(data()$edges$testedByDapSeq)/dim(data()$edges)[1],3)*100,
            "% (liens rouges + noirs). Parmis elles, ", round(sum(data()$edges$DapSeqAproved) / sum(data()$edges$testedByDapSeq),3)*100 ,
            "% ont été démontrées par DAP-Seq (liens rouges).")
    })
    
    
    observeEvent(input$nlpTargets, {
      
      if(input$nlpTargets){
        data <- data()
        data$nodes$group <- data$nodes$NLP_Target
        visNetworkProxy("network")%>% visUpdateNodes(data$nodes) %>%visGroups(groupname = "1", color = "darkred")
      }
      else{
        visNetworkProxy("network")  %>%  visUpdateNodes(data()$nodes)  %>% visGroups(groupname = "Regulator", 
                                                                                     size = 28, color = list("background" = "#003399", "border"="#FFFFCC"), 
                                                                                     shape = "square") %>% 
          visGroups(groupname = "Target Gene", color = "#77EEAA") 
      }
    })
}


# Run the application 
shinyApp(ui = ui, server = server)