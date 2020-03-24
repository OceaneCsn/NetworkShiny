library(PLNmodels)
library(ggplot2)
library(igraph)
library(GENIE3)
library(knitr)
PLN_network <- function(data, DEGenes, plot_path=F){
  # covariables
  groups <- str_split_fixed(colnames(data), "_", 2)[,1]
  co2 <- str_split_fixed(groups, "", 3)[,1]
  nitrate <- factor(str_split_fixed(groups, "", 3)[,2])
  nitrate <- relevel(nitrate, "N")
  fer <- factor(str_split_fixed(groups, "", 3)[,3])
  fer = relevel(fer, "F")
  covariates <- data.frame(row.names =colnames(data), co2,nitrate, fer)
  
  # preparation des données
  counts <- round(t(data[DEGenes,]), 0)
  plnData <- prepare_data(counts = counts, covariates = covariates)
  network_models <- PLNnetwork(Abundance ~ nitrate + fer + co2 +offset(log(Offset)), data = plnData)
  network_models
  network_models$criteria %>% head() %>% knitr::kable()
  plot(network_models, "diagnostic")
  plot(network_models)
  if(plot_path==T){
    coefficient_path(network_models, corr = TRUE) %>% 
      ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) + 
      geom_line(show.legend = FALSE) +  coord_trans(x="log10") + theme_bw()
  }
  
  model_StARS <- getBestModel(network_models, "StARS")
  print(model_StARS)
  net <- plot(model_StARS)
  plot.igraph(net)
  #plot.igraph(net, vertex.size = 10, vertex.label.cex = 0.5) 
  plot(model_StARS, type = "support", output = "corrplot")
  
  # Verification des predictions du modele
  data <- data.frame(
    fitted   = as.vector(fitted(model_StARS)),
    observed = as.vector(counts)
  ) 
  print(ggplot(data, aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme_bw() + annotation_logticks())
  return(net)
  
  netStats(g)
}

genie <- function(normalized.count, regressors, targets, nTrees=1000, nCores=5, top = 0.05, fixedLinkNumber = NA,
                  returnLinks = F){
  print(head(normalized.count))
  print(intersect(rownames(normalized.count),regressors))
  print(targets)
  mat <- GENIE3(normalized.count, regulators = intersect(rownames(normalized.count),regressors), targets = targets, treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = nCores,verbose = T)
  #hist(mat)
  
  return(mat)
}

getNet <- function(mat, fixedLinkNumber=NA, top=0.05){
  links <- getLinkList(mat)
  if(is.na(fixedLinkNumber)) {links <- links[1:round(dim(links)[1]*top,0),]}
  else {links <- links[1:fixedLinkNumber,]}
  print(paste0("Number of links : ", dim(links)[1]))
  g <- graph.data.frame(links, directed = F)
  return(g)
}

networkData <- function(net, ontologies){
  ont <- ontologies[match(V(net)$name, ontologies$ensembl_gene_id),]
  data <- toVisNetworkData(net)
  
  #attributs des noeuds et des liens
  data$nodes$Ontology <- ont[match(data$nodes$id, ont$ensembl_gene_id),]$external_gene_name
  data$nodes$description <-ont[match(data$nodes$id, ont$ensembl_gene_id),]$description
  data$nodes$group <- ifelse(data$nodes$id %in% TF$AGI, "Regulator", "Target Gene")
  
  importance <- (degree(net)/max(degree(net))+betweenness(net)/max(betweenness(net)))*0.5
  
  data$nodes$ranking <- importance[match(data$nodes$id, names(importance))]
  data$edges$value <- data$edges$weight
  #plot
  
  print(kable(head(data$nodes[order(-data$nodes$ranking),], n=40)))
  
  return(data)
}

plotNetwork <- function(data){
  visNetwork(nodes = data$nodes, edges = data$edges)%>% 
    #visEdges(smooth = FALSE, color = '#333366') %>% 
    visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
    visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10, 
               maxVelocity = 10, stabilization = F)%>%
    visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = TRUE)%>% 
    visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
    visGroups(groupname = "Regulator", size = 28,
              color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
    visGroups(groupname = "Target Gene", color = "#77EEAA") %>% 
    visNodes(borderWidth=0.5) 
}

# getName <- function(gene, ontologies){
#   if(is.na)
# }

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
