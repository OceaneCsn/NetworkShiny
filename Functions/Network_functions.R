library(ggplot2)
library(igraph)
library(knitr)
library(visNetwork)

plotNetwork <- function(data){
  visNetwork(nodes = data$nodes, edges = data$edges)%>% 
    visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
    visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10, 
               maxVelocity = 10, stabilization = F)%>%
    visOptions(selectedBy = "group", highlightNearest = F, nodesIdSelection  = TRUE, collapse = F)%>% 
    visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
    visGroups(groupname = "Regulator", size = 28,
              color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
    visGroups(groupname = "Target Gene", color = "#77EEAA") %>% 
    visNodes(borderWidth=0.5, font=list("size"=36)) 
}

NitrateGenes <- function(genesK, nGenes, ontologies){
  res <- data.frame(Gene = genesK)
  for(paper in colnames(nGenes)){
    res[,paper] <- ifelse(res$Gene %in% toupper(nGenes[,paper]), 1, 0)
  }
  res[,c("Name", "Description")] <- ontologies[match(res$Gene, ontologies$ensembl_gene_id),c("external_gene_name", "description")]
  res$NitrateScore <- rowSums(res[,grepl("_", colnames(res))])
  res <- res[order(-res$NitrateScore),]
  return(res)
}

networkData <- function(net, ontologies, TF){
  ont <- ontologies[match(V(net)$name, ontologies$ensembl_gene_id),]
  data <- toVisNetworkData(net)
  
  #attributs des noeuds et des liens
  data$nodes$Ontology <- ont[match(data$nodes$id, ont$ensembl_gene_id),]$external_gene_name
  data$nodes$description <-ont[match(data$nodes$id, ont$ensembl_gene_id),]$description
  data$nodes$group <- ifelse(data$nodes$id %in% TF$AGI, "Regulator", "Target Gene")
  
  importance <- (degree(net)/max(degree(net))+betweenness(net)/max(betweenness(net)))*0.5
  
  data$nodes$ranking <- importance[match(data$nodes$id, names(importance))]
  data$edges$value <- data$edges$weight
  print(kable(head(data$nodes[order(-data$nodes$ranking),], n=40)))
  
  return(data)
}


getTargets <- function(tf, data){
  targets <- data$edges[data$edges$from==tf,"to"]
  names <- data$nodes[data$nodes$id %in% targets, "label"]
  return(paste(names, collapse = ', '))
}

getRegulators <- function(gene, data){
  regs <- data$edges[data$edges$to==gene,"from"]
  names <- data$nodes[data$nodes$id %in% regs, "label"]
  return(paste(names, collapse = ', '))
}

getTargetsNumber <- function(tf, data){
  return(length(data$edges[data$edges$from==tf,"to"]))
}

colorEdge <- function(tested, validated) {
  if (tested & validated)
    return("red")
  if (tested & !validated)
    return("black")
  if (!tested & !validated)
    return("white")
}