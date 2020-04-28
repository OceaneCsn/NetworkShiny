# personnal functions
source("./Functions/Visu.R")
source("./Functions/Network_functions.R")
source("./Functions/ClusteringFunctions.R")

# expression levels and ontologies
load("./Data/OntologyAllGenes.RData")
load("./Data/normalized.count_At.RData")
load("./Data/NitrateGenes.RData")
load("./Data/DEGsListsFiltered.RData")
load("./Data/PlnTFBDRegulatorsList.RData")
load("./Data/N_uptake_genes.RData")

load(paste0("./NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData"))
load("./NetworkReference/Gaudinier_Nature_2018_TFs_Nitrates_Ref.RData")
#load("./NetworkReference/Gaudinier_Correlations_Litterature.RData")

networkCoseqMatching <- list("CO2DEGenes_faibleNitrate_CO2-N.RData" = "AmbientCO2_LowNitrateFe-ElevatedCO2_LowNitrateFeNoIronStarv.RData",
                             "CO2DEGenes_IronStarv_CO2-Fe.RData" = "AmbientCO2_HighNitrate_FeStarvation-ElevatedCO2_HighNitrate_FeStarvation_NoNitrateStarv.RData",
                             "CO2DEGenes_IronStarv_LowNitrate_CO2-N-Fe.RData" = "AmbientCO2_LowNitrate_FeStarvation-ElevatedCO2_LowNitrate_FeStarvation.RData")


universe <- ontologies$entrezgene_id

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
  data <- reactive({
    load(paste0("./NetworkData/", input$select))
    data$nodes$label <- data$nodes$Ontology
    data$edges$color <- '#333366'
    data$edges$Regulator_Name <-
      ontologies[match(data$edges$from, ontologies$ensembl_gene_id), ]$external_gene_name
    data$nodes$description <-
      ifelse(
        grepl("\\[", data$nodes$description),
        str_split_fixed(data$nodes$description, "\\[", 2)[, 1],
        data$nodes$description
      )
    data$edges$Target_Name <-
      ontologies[match(data$edges$to, ontologies$ensembl_gene_id), ]$external_gene_name
    data
  })
  
  
  
  tfs <- reactive({
    return(data()$nodes[data()$nodes$group == "Regulator", "id"])
  })
  
  targets <- reactive({
    return(data()$nodes[data()$nodes$group == "Target Gene", "id"])
  })
  
  
  output$tfNumber <- renderValueBox({
    valueBox(value = length(tfs()), subtitle = " Regulators",
             color = "navy")
  })
  
  
  output$targetNumber <- renderValueBox({
    valueBox(value = length(targets()), subtitle = "Target Genes",
    color = "teal")
  })
  
  
  
  output$network <- renderVisNetwork({
    visNetwork(nodes = data()$nodes, edges = data()$edges) %>%
      visEdges(smooth = FALSE,
               arrows = 'to',
               color = '#333366') %>%
      # visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10,
      #            maxVelocity = 10, stabilization = F)%>%
      visPhysics(
        solver = "forceAtlas2Based",
        timestep = 0.6,
        minVelocity = 12,
        maxVelocity = 10,
        stabilization = F
      ) %>%
      visOptions(
        selectedBy = "group",
        highlightNearest = TRUE,
        nodesIdSelection  = TRUE,
        collapse = F
      ) %>%
      visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>%
      visGroups(
        groupname = "Regulator",
        size = 28,
        color = list("background" = "#003399", "border" = "#FFFFCC"),
        shape = "square"
      ) %>%
      visGroups(groupname = "Target Gene",
                color = list("background" = "#77EEAA", hover = "grey")) %>%
      visNodes(borderWidth = 0.5, font = list("size" = 35))  %>% visInteraction(hover = TRUE)
    
  })
  
  observe({
    visNetworkProxy("network") %>% visFocus(id = input$geneVis, scale = 1) %>% visSelectNodes(id = input$geneVis)
  })
  
  output$SelectedGene <- renderPrint({
    print(input$click)
  })
  
  output$Ontologies <- DT::renderDataTable({
    if (is.null(input$click)) {
      data()$nodes
    }
    else{
      neighboors <-
        unique(union(data()$edges[grepl(input$click[1], data()$edges$from), ]$to,
                     data()$edges[grepl(input$click[1], data()$edges$to), ]$from))
      data()$nodes[c(input$click[1], neighboors), ]
    }
  })
  
  output$nitrate <- DT::renderDataTable({
    if (is.null(input$click)) {
      res <- NitrateGenes(data()$nodes$id, nGenes, ontologies)
      res[, !grepl("_", colnames(res))]
    }
    else{
      neighboors <-
        unique(union(data()$edges[grepl(input$click[1], data()$edges$from), ]$to,
                     data()$edges[grepl(input$click[1], data()$edges$to), ]$from))
      res <-
        NitrateGenes(c(input$click[1], neighboors), nGenes, ontologies)
      res[, !grepl("_", colnames(res))]
    }
  })
  
  output$Ranking <- DT::renderDataTable({
    net <- igraph::graph_from_data_frame(d = data()$edges)
    importance <-
      (degree(net) / max(degree(net)) + betweenness(net) / max(betweenness(net))) *
      0.5
    
    data <- data()
    
    data$nodes$Betweenness <-
      betweenness(net)[match(data$nodes$id, names(betweenness(net)))]
    data$nodes$Degree <-
      degree(net)[match(data$nodes$id, names(degree(net)))]
    data$nodes$Centrality <-
      importance[match(data$nodes$id, names(importance))]
    data$nodes$Rank <- order(-data$nodes$Centrality)
    data$nodes[order(-data$nodes$Centrality), c("description",
                                                "Ontology",
                                                "Centrality",
                                                "Betweenness",
                                                "Degree")]
    
  })
  
  output$expression_plot <- renderPlot({
    getExpression(input$click, normalized.count = normalized.count)
    
  })
  
  
  output$NetStats <- renderPlot({
    net <- igraph::graph_from_data_frame(d = data()$edges)
    netStats(net, TFs = tfs(), targets = targets())
  })
  
  
  observeEvent(input$colorFromNitrate, {
    newNodes <- data()$nodes
    nScores <- NitrateGenes(data()$nodes$id, nGenes, ontologies)
    newNodes$group <-
      nScores[match(newNodes$id, nScores$Gene), "NitrateScore"]
    visNetworkProxy("network") %>%
      visUpdateNodes(nodes = newNodes) %>% visGroups(groupname = "6",
                                                     size = 28,
                                                     color = "#660000") %>%
      visGroups(groupname = "5", color = "#990000") %>%
      visGroups(groupname = "4", color = "#cc0000") %>%
      visGroups(groupname = "3", color = "#e06666") %>%
      visGroups(groupname = "2", color = "#ea9999") %>%
      visGroups(groupname = "1", color = "#f4cccc") %>%
      visGroups(groupname = "0", color = "#ffffff") %>%
      visGroups(groupname = "7", color = "#330000")
  })
  
  
  output$heatmap <- renderPlotly({
    if (is.null(input$click)) {
      heatmapPerso(
        normalized.count,
        conds = "all",
        genes = c("AT1G01150", "AT1G01590")
      )
    }
    else{
      neighboors <-
        unique(union(data()$edges[grepl(input$click[1], data()$edges$from), ]$to,
                     data()$edges[grepl(input$click[1], data()$edges$to), ]$from))
      heatmapPerso(
        normalized.count,
        conds = "all",
        genes = c(input$click[1], neighboors)
      )
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
    for (comp in names(DEGs)) {
      if (input$geneToSearch %in% DEGs[[comp]])
        print(comp)
    }
    
  })
  
  output$commonLinks <- DT::renderDataTable({
    load(paste0("./NetworkData/", input$select))
    data$edges$pairs <- paste(data$edges$from, data$edges$to)
    gaudinier$edges$pairs <-
      paste(gaudinier$edges$from, gaudinier$edges$to)
    
    commonLinks <-
      intersect(gaudinier$edges$pairs, data$edges$pairs)
    data$edges[data$edges$pairs %in% commonLinks, c("from", "to")]
  })
  
  output$commonNodes <- DT::renderDataTable({
    load(paste0("./NetworkData/", input$select))
    data$nodes[data$nodes$id %in% gaudinier$nodes$id, c("Ontology", "description", "group", "ranking")]
  })
  
  output$netIntersection <- renderVisNetwork({
    load(paste0("./NetworkData/", input$select1))
    network1 <- igraph::graph_from_data_frame(d = data$edges)
    
    load(paste0("./NetworkData/", input$select2))
    network2 <- igraph::graph_from_data_frame(d = data$edges)
    
    inter <-
      igraph::intersection(network1, network2, keep.all.vertices = F)
    dataInter <- networkData(inter, ontologies, TF)
    dataInter$nodes$label <- dataInter$nodes$Ontology
    
    nScores <-
      NitrateGenes(dataInter$nodes$id, nGenes, ontologies)
    dataInter$nodes$group <-
      nScores[match(dataInter$nodes$id, nScores$Gene), "NitrateScore"]
    
    
    visNetwork(nodes = dataInter$nodes, edges = dataInter$edges) %>%
      visEdges(smooth = FALSE,
               arrows = 'to',
               color = '#333366') %>%
      visPhysics(
        solver = "forceAtlas2Based",
        timestep = 0.6,
        minVelocity = 12,
        maxVelocity = 10,
        stabilization = F
      ) %>%
      visOptions(
        selectedBy = "group",
        highlightNearest = TRUE,
        nodesIdSelection  = TRUE,
        collapse = F
      ) %>%
      visEvents(click = "function(nodes){
                Shiny.onInputChange('click', nodes.nodes);
                ;}") %>%  visNodes(borderWidth = 0.5, font = list("size" =
                                                                    35))  %>% visInteraction(hover = TRUE) %>% visGroups(groupname = "6",
                                                                                                                         size = 28,
                                                                                                                         color = "#660000") %>%
      visGroups(groupname = "5", color = "#990000") %>%
      visGroups(groupname = "4", color = "#cc0000") %>%
      visGroups(groupname = "3", color = "#e06666") %>%
      visGroups(groupname = "2", color = "#ea9999") %>%
      visGroups(groupname = "1", color = "#f4cccc") %>%
      visGroups(groupname = "0", color = "#ffffff") %>%
      visGroups(groupname = "7", color = "#330000")
    
    
  })
  
  output$commonNodes12 <- DT::renderDataTable({
    load(paste0("./NetworkData/", input$select1))
    data1 <- data
    
    load(paste0("./NetworkData/", input$select2))
    
    data1$nodes[data1$nodes$id %in% data$nodes$id, c("Ontology", "description", "group", "ranking")]
  })
  
  output$commonLinks12 <- DT::renderDataTable({
    load(paste0("./NetworkData/", input$select1))
    network1 <- igraph::graph_from_data_frame(d = data$edges)
    
    load(paste0("./NetworkData/", input$select2))
    network2 <- igraph::graph_from_data_frame(d = data$edges)
    
    inter <-
      igraph::intersection(network1, network2, keep.all.vertices = F)
    dataInter <- networkData(inter, ontologies, TF)
    
    dataInter$edges[, c("from", "to")]
  })
  
  #################### Clustering network
  
  
  output$netClustering <- renderVisNetwork({
    plotNetwork(dataClust()) %>% visGroups(groupname = "-1", color = "lightgrey")
  })
  
  
  dataClust <- reactive({
    load(paste0("./NetworkData/", input$selectClusterNetwork))
    data$nodes$role <- data$nodes$group
    
    data$nodes$group <-
      data$nodes[, paste0(input$clustType, "Cluster")]
    updateSelectInput(session, "Module", choices = data$nodes$group)
    data$nodes$label <- data$nodes$Ontology
    data$edges$color <- '#333366'
    data$nodes$description <-
      ifelse(
        grepl("\\[", data$nodes$description),
        str_split_fixed(data$nodes$description, "\\[", 2)[, 1],
        data$nodes$description
      )
    data
  })
  
  observe({
    visNetworkProxy("netClustering") %>% visOptions(selectedBy = list(variable = "group", selected = input$Module))
  })
  
  output$communityList <- DT::renderDataTable({
    dataClust()$nodes[dataClust()$nodes$group == input$Module, ]
  })
  
  output$GOEnrich <- renderPlotly({
    ids <-
      as.character(ontologies[match(dataClust()$nodes[dataClust()$nodes$group ==
                                                        input$Module, ]$id, ontologies$ensembl_gene_id), ]$entrezgene_id)
    withProgress(message = "Ontologies Enrichment", {
      simOnt <- OntologyEnrich(ids, as.character(universe))
    })
    simOnt@result <-
      simOnt@result[order(-simOnt@result$p.adjust), ]
    values <- str_split_fixed(simOnt@result$GeneRatio, "/", 2)
    simOnt@result$GeneRatio <-
      as.numeric(values[, 1]) / as.numeric(values[, 2])
    ggplotly(
      ggplot(
        data = simOnt@result,
        aes(x = Description, y = GeneRatio, fill = p.adjust)
      ) + geom_bar(stat = "identity") + coord_flip() +
        ggtitle(paste(
          "Enriched Ontologies for module", input$Module
        )) + ylab("") + xlab("Gene Ratio") +
        theme(
          plot.title = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 15, angle = 20),
          axis.text.x = element_text(
            size = 8,
            angle = 0,
            hjust = 0,
            colour = "grey50"
          ),
          axis.title.y = element_blank()
        )
    )
  })
  
  output$GOEnrichComp <- renderPlot({
    idsList <- list()
    for (k in unique(dataClust()$nodes$group)) {
      idsList[[as.character(k)]] <-
        na.omit(as.character(ontologies[match(dataClust()$nodes[dataClust()$nodes$group ==
                                                                  k, ]$id, ontologies$ensembl_gene_id), ]$entrezgene_id))
    }
    withProgress(message = 'Ontologies enrichment comparison', {
      compareOnt(idsList = idsList, as.character(universe))
    })
  })
  
  output$topTFs <- DT::renderDataTable({
    dataK <-
      dataClust()$nodes[dataClust()$nodes[, paste0(input$clustType, "Cluster")] ==
                          input$Module &
                          dataClust()$nodes$role == "Regulator", ]
    top <- dataK[order(-dataK$ranking), ]
    top <- top[1:round(dim(dataK)[1] * input$topRate / 100, 0), ]
    top$targetNumber <-
      sapply(top$id, getTargetsNumber, dataClust())
    top$targets <- sapply(top$id, getTargets,
                          dataClust())
    top[, c("label",
            "description",
            "ranking",
            "targetNumber",
            "targets",
            "NitrateScore")]
  })
  
  cluster <- reactive({
    load(paste0("./ClusteringData/", networkCoseqMatching[[input$selectClusterNetwork]]))
    cluster
  })
  
  output$profiles <- renderPlot({
    
    
    if(input$allClustersProfile){plotProfileFromNetwork(dataClust(), normalized.count,
                                                        paste0(input$clustType, "Cluster"), boxplot = input$boxplot,
                                                        removeIronStarv = grepl("NoIronStarv", networkCoseqMatching[[input$selectClusterNetwork]]),
                                                        removeNitrateStarv = grepl("NoNitrateStarv", networkCoseqMatching[[input$selectClusterNetwork]]))}
    else{plotProfileFromNetwork(dataClust(), normalized.count,
                                paste0(input$clustType, "Cluster"), k = input$Module, boxplot=input$boxplot,
                                removeIronStarv = grepl("NoIronStarv", networkCoseqMatching[[input$selectClusterNetwork]]),
                                removeNitrateStarv = grepl("NoNitrateStarv", networkCoseqMatching[[input$selectClusterNetwork]]))}
  })
  
  ########## glm functions 
  
  glm <- reactive({
  if(input$clustType=="coseq"){
    glmCluster(DEgenes = names(cluster()[[1]][cluster()[[1]]==input$Module]), 
               normalized.count = data.frame(cluster()[[2]]@tcounts)) 
  }
  else{
    DEgenes = dataClust()$nodes[dataClust()$nodes[,paste0(input$clustType, "Cluster")]==input$Module, "id"]
    glmCluster(DEgenes, 
               normalized.count = normalized.count, removeIronStarv = grepl("NoIronStarv", networkCoseqMatching[[input$selectClusterNetwork]]),
               removeNitrateStarv = grepl("NoNitrateStarv", networkCoseqMatching[[input$selectClusterNetwork]]))
  }
  
  })
  
  output$coefsPlot <- renderPlotly({
    print(glm())
    plotGlmCluster(glm())
  })
  
  output$summary <- renderPrint({
    summary(glm())
  })
  
  ############################################ fav Genes
  
  geneList <- reactive({
    genes <- as.vector(input$groupGenes)
    res <- c()
    for (gene in genes) {
      res <-
        c(res, dataClust()$nodes$label[grepl(gene, dataClust()$nodes$label)])
    }
    res
  })
  
  
  output$favGenes <- DT::renderDataTable({
    top <- dataClust()$nodes[dataClust()$nodes$label %in% geneList(), ]
    
    top$targets <- sapply(top$id, getTargets,
                          dataClust())
    
    top$regulators <- sapply(top$id, getRegulators,
                             dataClust())
    
    top <- top[order(top[,paste0(input$clustType, "Cluster")]), c(
      "label",
      "description",
      "regulators",
      "targets",
      "louvainCluster",
      "coseqCluster",
      "consensusCluster"
    )]
    formatStyle(datatable(top), columns = paste0(c("coseq", "louvain", "consensus"), "Cluster"), 
                       target = c("cell", "row"),
                       backgroundColor = styleEqual(c(1:12), c("#ccccff","#99ccff", "#ffcc00",
                                                    "#ffffcc", "#ccffff","#ccffcc",
                                                    "#ffe6e6", "#ffd9b3", "#ffcce6",
                                                    "#f2ccff", "#33adff", " #00cc00")))
  })
  
  output$rootGenes <- DT::renderDataTable({
    
    geneList <- c("AT5G65230", "AT5G16770", "AT1G34670", "AT5G16770", "AT1G58100",
                  "AT4G11880", "AT4G38620", "AT3G19580", "AT3G03660")
    top <- dataClust()$nodes[dataClust()$nodes$id %in% geneList, ]
    
    top$targets <- sapply(top$id, getTargets,
                          dataClust())
    
    top$regulators <- sapply(top$id, getRegulators,
                             dataClust())
    
    top <- top[order(top[,paste0(input$clustType, "Cluster")]), c(
      "label",
      "description",
      "regulators",
      "louvainCluster",
      "coseqCluster",
      "consensusCluster"
    )]
    formatStyle(datatable(top), columns = paste0(c("coseq", "louvain", "consensus"), "Cluster"), 
                target = c("cell", "row"),
                backgroundColor = styleEqual(c(1:12), c("#ccccff","#99ccff", "#ffcc00",
                                                        "#ffffcc", "#ccffff","#ccffcc",
                                                        "#ffe6e6", "#ffd9b3", "#ffcce6",
                                                        "#f2ccff", "#33adff", " #00cc00")))
  })
  
    output$nitrateUptakeGenes <- DT::renderDataTable({
      
      top <- dataClust()$nodes[dataClust()$nodes$id %in% n_uptake_genes, ]
      
      top$targets <- sapply(top$id, getTargets,
                            dataClust())
      
      top$regulators <- sapply(top$id, getRegulators,
                               dataClust())
      
      top <- top[order(top[,paste0(input$clustType, "Cluster")]), c(
        "label",
        "description",
        "regulators",
        "targets",
        "louvainCluster",
        "coseqCluster",
        "consensusCluster"
      )]
      formatStyle(datatable(top), columns = paste0(c("coseq", "louvain", "consensus"), "Cluster"), 
                  target = c("cell", "row"),
                  backgroundColor = styleEqual(c(1:12), c("#ccccff","#99ccff", "#ffcc00",
                                                          "#ffffcc", "#ccffff","#ccffcc",
                                                          "#ffe6e6", "#ffd9b3", "#ffcce6",
                                                          "#f2ccff", "#33adff", " #00cc00")))
    
  })
  
  #################################### DAPSeq functions
  
  observeEvent(input$edgesDap, {
    if (input$edgesDap) {
      data <- data()
      data$edges$color <-
        mapply(colorEdge,
               data$edges$testedByDapSeq,
               data$edges$DapSeqAproved)
      visNetworkProxy("network") %>% visUpdateEdges(data$edges) %>%
        visGroups(
          groupname = "Regulator",
          size = 28,
          color = list("background" = "grey", "border" = "#CCCCCC"),
          shape = "square"
        ) %>%
        visGroups(groupname = "Target Gene", color = "lightgrey") %>%
        visNodes(borderWidth = 0.5, font = list("size" = 36))
      
    }
    
    else{
      visNetworkProxy("network")  %>%  visUpdateEdges(data()$edges) %>% visGroups(
        groupname = "Regulator",
        size = 28,
        color = list("background" = "#003399", "border" =
                       "#FFFFCC"),
        shape = "square"
      ) %>%
        visGroups(groupname = "Target Gene", color = "#77EEAA")
    }
    
  })
  
  
  
  output$edgesListDap <- DT::renderDataTable({
    dataDap <- data()$edges[data()$edges$DapSeqAproved, ]
    if (is.null(input$click)) {
      dataDap[, c("from", "Regulator_Name", "to", "Target_Name")]
    }
    else{
      dataDap[grepl(input$click, dataDap$from), c("from", "Regulator_Name", "to", "Target_Name")]
    }
    
  })
  
  output$dapStats <- renderText({
    paste(
      "Interactions comportant un TF etudie par DAP-Seq : ",
      round(sum(data()$edges$testedByDapSeq) / dim(data()$edges)[1], 3) * 100,
      "% (liens rouges + noirs). Parmis elles, ",
      round(
        sum(data()$edges$DapSeqAproved) / sum(data()$edges$testedByDapSeq),
        3
      ) * 100 ,
      "% ont été démontrées par DAP-Seq (liens rouges)."
    )
  })
  
  
  observeEvent(input$nlpTargets, {
    if (input$nlpTargets) {
      data <- data()
      data$nodes$group <- data$nodes$NLP_Target
      visNetworkProxy("network") %>% visUpdateNodes(data$nodes) %>% visGroups(groupname = "1", color = "darkred")
    }
    else{
      visNetworkProxy("network")  %>%  visUpdateNodes(data()$nodes)  %>% visGroups(
        groupname = "Regulator",
        size = 28,
        color = list("background" = "#003399", "border" = "#FFFFCC"),
        shape = "square"
      ) %>%
        visGroups(groupname = "Target Gene", color = "#77EEAA")
    }
  })
}


# Run the application
#shinyApp(ui = ui, server = server)