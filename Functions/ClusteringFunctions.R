library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
library(gridExtra)
library(reshape2)
suppressMessages(library(coseq, warn.conflicts = F, quietly = T))
library(plotly)
library(stringr)


####################### Ontology analysis


OntologyEnrich <- function(ids, universe, plot = T, simCutoff = 0.8){
  # ids and universe must be entrez gene ids
  ego <- enrichGO(gene = ids,
                  OrgDb = org.At.tair.db,
                  ont = "BP",
                  universe = universe,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable = TRUE)
  # Elimine les redondances, en fusionnant les GO terms dont la similarite excede le cutoff
  simpOnt <- clusterProfiler::simplify(ego, cutoff=simCutoff, by="p.adjust", select_fun=min)
  #print(barplot(simpOnt, showCategory = 40, font.size = 10))
  #print(emapplot(simpOnt, font.size = 30, layout = "kk"))
  return(simpOnt)
}

compareOnt <- function(idsList, universe, simCutoff = 0.8){
  ckreg <- compareCluster(geneCluster = idsList, fun = "enrichGO", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", 
                          pvalueCutoff = 0.01, qvalueCutoff = 0.05, universe = as.character(universe))
  ckreg@compareClusterResult
  simCk <- clusterProfiler::simplify(ckreg, cutoff=simCutoff, by="p.adjust", select_fun=min)
  resCk <- simCk@compareClusterResult
  print(clusterProfiler::dotplot(simCk, x = ~Cluster, showCategory = 30, font.size = 15))
}


################ profiles


translate <- function(text){
  res = ""
  if(grepl("c", text)){res = paste0(res, "AmbientCO2,")}
  else{res = paste0(res, "ElevatedCO2,")}
  if(grepl("N", text)){res = paste0(res, "HightNitrate,")}
  else{res = paste0(res, "LowNitrate,")}
  if(grepl("f", text)){res = paste0(res, "IronStarvation")}
  else{res = paste0(res, "IronSupply")}
  return(res)
}

plotProfile <- function(cluster, k="none", boxplot=T, expression = "profiles"){
  # plot all the profiles or the profile of cluster k
  results <- cluster[[2]]
  clusters <- cluster[[1]]
  if(expression=="profiles"){profiles <- data.frame(results@y_profiles)
  ylab <- "Normalized expression/Mean(Normalized expression)"}
  if(expression=="counts") {
    profiles <- data.frame(log(as.matrix(results@tcounts)+1))
    ylab <- "log(Normalized expresion)"
  }
  profiles$gene <- rownames(profiles)
  d <- melt(profiles)
  d$group <- str_split_fixed(d$variable, '_', 2)[,1]
  d$cluster <- clusters[match(d$gene, names(clusters))]
  d$geneRep <- paste0(d$gene, substr(d$variable,4,5))
  if(k=="none"){
    g <- ggplot(data = d, aes(x=group, y=value))  +facet_wrap(~cluster, nrow=3) 
  }
  else{
    g <- ggplot(data = d[d$cluster==k,], aes(x=group, y=value))
  }
  if(boxplot) g <- g + geom_boxplot(alpha=0.7, lwd=1.2, aes( color = group), fill = "grey", outlier.color = "black",outlier.alpha =0.1)  + geom_jitter(width = 0.1, alpha=0.0015)
  else{
    g <- g+ geom_line(alpha=0.09,lwd=1.2, color="#333366", aes(group=geneRep))
  }
  
  g <- g +theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 20),legend.position="bottom",
                legend.title = element_text(size = 2, face="bold"), legend.text = element_text(size=15, angle=0),
                axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 0, hjust = 0, colour = "grey50"),legend.text.align=1,
                axis.title=element_text(size=24)) + xlab("") + ylab(ylab) + scale_colour_discrete("", labels=sapply(levels(as.factor(d$group)),translate)) +
    stat_summary(fun.y=median, geom="line", aes(group=1), alpha=0.1, size = 1.5) 
  if(expression=="profiles") g <- g + ylim(0, 0.25) 
  g
}

plotProfileFromNetwork <- function(netData, normalized.count, clustType, k="none", boxplot=T, expression = "profiles",
                                   removeIronStarv = F, removeNitrateStarv=F){
  # plot all the profiles or the profile of cluster k
  if(removeIronStarv) normalized.count <- normalized.count[,grepl("F", colnames(normalized.count))]
  if(removeNitrateStarv) normalized.count <- normalized.count[,grepl("N", colnames(normalized.count))]
  
  normalized.count <- normalized.count[netData$nodes$id,]
  if(expression=="profiles"){profiles <- data.frame(normalized.count/rowSums(normalized.count))
    ylab <- "Normalized expression/Mean(Normalized expression)"}
  if(expression=="counts") {
    profiles <- data.frame(log(as.matrix(normalized.count)+1))
    ylab <- "log(Normalized expresion)"
  }
  profiles$gene <- rownames(profiles)
  d <- melt(profiles)
  print(head(d))
  d$group <- str_split_fixed(d$variable, '_', 2)[,1]
  d$cluster <- netData$nodes[match(d$gene, netData$nodes$id),clustType]
  d$geneRep <- paste0(d$gene, substr(d$variable,4,5))
  print(head(d))
  if(k=="none"){
    g <- ggplot(data = d, aes(x=group, y=value))  +facet_wrap(~cluster, nrow=3) 
  }
  else{
    g <- ggplot(data = d[d$cluster==k,], aes(x=group, y=value))
  }
  if(boxplot) g <- g + geom_boxplot(alpha=0.7, lwd=1.2, aes( color = group), fill = "grey", outlier.color = "black",outlier.alpha =0.1)  + geom_jitter(width = 0.1, alpha=0.0015)
  else{
    g <- g+ geom_line(alpha=0.09,lwd=1.2, color="#333366", aes(group=geneRep))
  }
  
  g <- g +theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 20),legend.position="bottom",
                legend.title = element_text(size = 2, face="bold"), legend.text = element_text(size=15, angle=0),
                axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 0, hjust = 0, colour = "grey50"),legend.text.align=1,
                axis.title=element_text(size=24)) + xlab("") + ylab(ylab) + scale_colour_discrete("", labels=sapply(levels(as.factor(d$group)),translate)) +
    stat_summary(fun.y=median, geom="line", aes(group=1), alpha=0.1, size = 1.5) 
  if(expression=="profiles") g <- g + ylim(0, 0.25) 
  g
}



###################### Poisson GLM

glmCluster <- function(DEgenes, normalized.count, removeIronStarv = F,
                       removeNitrateStarv=F){
  if(removeIronStarv) normalized.count <- normalized.count[,grepl("F", colnames(normalized.count))]
  if(removeNitrateStarv) normalized.count <- normalized.count[,grepl("N", colnames(normalized.count))]
  
  glmData <- melt(round(normalized.count[DEgenes,], 0))
  print(head(glmData))
  if(grepl("AT", glmData[1,1])){glmData <- glmData[,2:length(colnames(glmData))]}
  colnames(glmData) <- c("Condition", "Counts")
  glmData <- glmData[sample(rownames(glmData)),]
  
  groups <- str_split_fixed(glmData$Condition, "_", 2)[,1]
  glmData$Co2 <- str_split_fixed(groups, "", 3)[,1]
  glmData$nitrate <- str_split_fixed(groups, "", 3)[,2]
  glmData$fer <- str_split_fixed(groups, "", 3)[,3]
  
  
  glmData$Co2 <- as.factor(ifelse(glmData$Co2 == "c", 0, 1))
  glmData$nitrate <- as.factor(ifelse(glmData$nitrate == "N", 0, 1))
  glmData$fer <- as.factor(ifelse(glmData$fer == "F", 0, 1))
  #glmData <- glmData[c("Counts", "Co2", "nitrate", "fer", "gene")]
  
  print(head(glmData))
  formula = "Counts ~ "
  
  for (factor in c("Co2", "nitrate", "fer")){
    if(length(levels(glmData[,factor]))>1){formula <- paste(formula, factor, '*')}
  }
  
  formula <- substr(formula, 1, nchar(formula)-1)
  
  glm <- glm(formula , data = glmData, family = poisson(link="log"))
  return(glm)
  
}

plotGlmCluster <- function(glm){
  coefs <- glm$coefficients[2:length(glm$coefficients)]
  d <- data.frame(Coefficient =str_remove_all(as.character(names(coefs)), "1"), Values = coefs)
  ggplotly(ggplot(d, aes(x=Coefficient, y=Values, fill=Coefficient)) + geom_bar(color = 'black', stat = "identity", alpha = 0.4) +
             theme(axis.text.x = element_text(size=25, angle = 320), legend.position = "none", 
                   axis.title.x = element_blank(), axis.title.y = element_blank(), 
                   plot.title = element_text(size=22, face="bold"), axis.text.y = element_text(size=20))+
             ggtitle("Generalized linear model's coefficients values"), tooltip = c("Values"))
}

