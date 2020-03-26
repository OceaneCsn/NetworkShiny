suppressMessages(library(stringr, warn.conflicts = F, quietly = T))
suppressMessages(library(tidyr, warn.conflicts = F, quietly = T))
suppressMessages(library(ggplot2, warn.conflicts = F, quietly = T))
suppressMessages(library(gridExtra, warn.conflicts = F, quietly = T))

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
#setwd("D:/These/CombinatoireRNASeqFeNCO2/")

heatmapPerso <- function(normalized.count, genes=NA, conds="all", specie="At", geneNames=NA, profiles=T){
  if(length(genes) < 1){genes <- rownames(normalized.count)[1:6]}
  #if(specie == "At") load("normalized.count_At.RData")
  #if(specie == "Sl") load("normalized.count_Sl.RData")
  if(profiles) normalized.count=normalized.count/rowSums(normalized.count)
  
  if (length(conds) ==1){
    conds = colnames(normalized.count)
  }else{conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))}
  df <- data.frame(t(normalized.count[genes,conds]))
  
  df$condition <- str_split_fixed(rownames(df), "_", 2)[,1]
  df$exactCondition <- rownames(df)
  data <- gather(data = df,key = gene,value = expression, -condition, -exactCondition)
  exp.heatmap <- ggplot(data = data, mapping = aes(x = exactCondition, y = gene,
                                                   fill = log(expression+0.1))) +
    geom_tile() + xlab(label = "Condition") +
    facet_grid(~ condition, switch = "x", scales = "free_x", space = "free_x") + labs(fill = "Normalized expression") +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank()) + scale_fill_distiller(palette = "YlGnBu") 
  if(length(geneNames) > 1) {
    exp.heatmap = exp.heatmap + scale_y_discrete(labels=geneNames)
  }
  return(ggplotly(exp.heatmap))
}



getExpression <- function(gene, normalized.count, conds = "all", specie = "At"){
  # Plots the expression levels of a given gene, using the normized.count data provoded.
  # conditions are all the columns of the data by default, or can be specified
  # biological replicated should be identified by _n

  if (length(conds) ==1){
    conds = colnames(normalized.count)
  }else{conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))}
  df <- normalized.count[gene, conds]
  library(reshape2)
  d<- melt(df, silent=T)
  d$group = str_split_fixed(rownames(d), "_", 2)[,1]
  
  p <- ggplot(data = d, aes(x=group, y=value, fill=group)) + geom_dotplot(binaxis = "y", stackdir = "center") +
    scale_fill_discrete(name = "Conditions", labels=sapply(levels(as.factor(d$group)),translate))+
    theme(strip.text.x = element_text(size = 26), plot.title = element_text(size=22, face="bold"),
          legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=20),
          axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 26, angle = 320, hjust = 0, colour = "grey50"),
          axis.title=element_text(size=17)) + ylab("Normalized counts") +
    ggtitle(paste("Normalized expression for ", gene)) + xlab("- C : elevated CO2 - c : ambiant CO2 - N : 10mM nitrate - n : 0.5mM nitrate - F : iron - f : iron starvation")
  return(p)
}

netStats <- function(g){
  degree<- degree(g)
  betweenness<- betweenness(g, weights=NA)
  Node_nw_st<- data.frame(degree, betweenness)
  deg <- ggplot( data = Node_nw_st, aes(x=degree)) +geom_histogram( bins = 50, fill="#69b3a2", color="#e9ecef", alpha=0.7) +
          ggtitle("Degree distribution") +
          theme(plot.title = element_text(size=15))
  bet <- ggplot( data = Node_nw_st, aes(x=betweenness)) +geom_histogram(bins = 50, fill="#E69F00", color="#e9ecef", alpha=0.7) +
          ggtitle("Betweeness distribution") +
          theme(plot.title = element_text(size=15))
  return(grid.arrange(deg, bet, nrow = 1))
}