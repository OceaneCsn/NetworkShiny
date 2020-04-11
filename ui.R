library(shiny)
library(shinythemes)
library(DT)
library(visNetwork)
library(stringr)
library(ggplot2)
library(shinydashboard)
library(igraph)
library(plotly)
if (version$os == "linux-gnu")
  options(bitmapType = 'cairo')


load("./Data/DEGsListsFiltered.RData")


listFiles <- list.files("./NetworkData/", full.names = F)
names(listFiles) = listFiles
files <- lapply(split(listFiles, names(listFiles)), unname)


ui <- dashboardPage(
  skin = "black",
  
  
  dashboardHeader(title = "Network visualisation"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "Network",
        tabName = "Network",
        icon = icon("project-diagram")
      ),
      menuItem(
        "Expression database",
        tabName = "Expression_data",
        icon = icon("seedling")
      ),
      menuItem(
        "Differential Expression database",
        tabName = "DEGs",
        icon = icon("table")
      ),
      menuItem(
        "Compare 2 networks",
        tabName = "Comparaison",
        icon = icon("greater-than-equal")
      ),
      menuItem(
        "Communities discovery",
        tabName = "Clustering",
        icon = icon("circle-notch")
      )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(
        ".shiny-progress {top: 50% !important;left: 50% !important;margin-top: -100px !important;margin-left: -250px !important; color: blue;font-size: 20px;font-style: italic;}"
      )
    ),
    tabItems(
      tabItem(
        tabName = "Network",
        selectInput(
          "select",
          label = h3("Select network data"),
          width = 700,
          choices = files,
          selected = "CO2DEGenes_faibleNitrate_CO2-N.RData"
        ),
        
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
            tabsetPanel(
              type = "tabs",
              tabPanel("Ontologies", DT::dataTableOutput("Ontologies")),
              tabPanel("Ontologies Nitrate", DT::dataTableOutput("nitrate")),
              tabPanel(
                "Heatmap",
                plotlyOutput("heatmap"),
                br(),
                plotOutput("expression_plot")
              ),
              tabPanel("Gene ranking", DT::dataTableOutput("Ranking")),
              tabPanel("Network statistics", plotOutput("NetStats")),
              tabPanel(
                "Comparison to gaudinier",
                h2(
                  "Combinatorial interactions between transcription factors and promoters of genes associated with nitrogen metabolism, signalling and nitrogen-associated processes"
                ),
                h3("Common links with state of the art network :"),
                DT::dataTableOutput("commonLinks"),
                h3('Common genes with state of the art network : '),
                DT::dataTableOutput("commonNodes")
              ),
              tabPanel(
                "DAP-Seq validated interactions",
                checkboxInput("edgesDap", "Color edges from DAP-Seq", FALSE),
                h3(
                  "To check if the predicted target's promoters have the regulator's binding motif"
                ),
                h2("Interactions validated by DAP-Seq :"),
                DT::dataTableOutput("edgesListDap", width = 700),
                h3("Taux de validation : "),
                textOutput("dapStats"),
                checkboxInput("nlpTargets", "Color NLP7 DAP-Seq targets", FALSE)
              )
            )
          )
        )
      ),
      tabItem(
        tabName = "Expression_data",
        textInput("gene", "Ask me a gene! (AGI)", value = "AT1G01020"),
        br(),
        plotOutput("expression_plot_specific", width = "1400px")
        
      ),
      
      tabItem(
        tabName = "DEGs",
        fixedRow(column(
          width = 6,
          h3(
            "See which genes aree differentially expressed in one transcriptome comparison :"
          ),
          selectInput(
            "comparison",
            label = h3("Select differentiel expression analysis"),
            choices = names(DEGs)
          ),
          DT::dataTableOutput("genesComp", width = "400px")
        )),
        column(
          width = 4,
          h3(
            "For a gene, find in which transcriptome comparison it is differentially expressed :"
          ),
          textInput("geneToSearch", "Ask me a gene! (AGI)", value = "AT1G08090"),
          verbatimTextOutput("comparisonsDEG")
        )
      ),
      
      tabItem(
        tabName = "Comparaison",
        
        selectInput(
          "select1",
          label = h3("Select first network data"),
          width = 600,
          choices = files,
          selected = "NitratesDEGenes_CO2-N.RData"
        ),
        
        selectInput(
          "select2",
          label = h3("Select second network data"),
          width = 600,
          choices = files,
          selected = "NitratesDEGenes_N.RData"
        ),
        
        fixedRow(
          column(
            width = 5,
            h3("Networks intersection :"),
            visNetworkOutput("netIntersection", height = "1000px")
          )
          ,
          column(
            width = 5,
            h2("Common links :"),
            DT::dataTableOutput("commonLinks12"),
            h2('Common genes : '),
            DT::dataTableOutput("commonNodes12")
          )
        )
      ),
      tabItem(
        tabName = "Clustering",
        
        column(
          width = 3,
          selectInput(
            "selectClusterNetwork",
            label = h3("Select network data"),
            width = 600,
            choices = files,
            selected = "CO2DEGenes_faibleNitrate_CO2-N.RData"
          )
        ),
        column(
          width = 3,
          selectInput(
            "clustType",
            label = h3("Select type of clustering"),
            width = 600,
            choices = list(
              "Topological clustering (Louvain)" = "louvain",
              "Expression clustering (COSEQ Poisson mixture models)" = "coseq"
            ),
            selected = "louvain"
          )
        ),
        
        column(
          width = 3,
          selectInput(
            "Module",
            label = h3("Select module"),
            width = 600,
            choices = NA
          )
        ),
        
        fixedRow(
          column(
            width = 4,
            h3("Network modules :"),
            visNetworkOutput("netClustering", height = "1000px")
          )
          ,
          column(
            width = 5,
            tabsetPanel(
              type = "tabs",
              tabPanel(
                "Genes list",
                h2("Genes in that Community :"),
                DT::dataTableOutput("communityList")
              ),
              tabPanel(
                "GO enrichment",
                h2("Ontologies enriched in this cluster : "),
                plotlyOutput("GOEnrich", width = 900, height = 900)
              ),
              tabPanel(
                "GO enrichment comparison",
                h2("Comparison of the communitites : "),
                plotOutput("GOEnrichComp", height = 900)
              ),
              tabPanel(
                "Regulators ranking",
                h2("Most influent regulators"),
                numericInput("topRate", label = "Top percentage of regulators to select : (%)", value = 20),
                DT::dataTableOutput("topTFs")
              )
            ),
            
          )
        )
      )
    )
  )
)