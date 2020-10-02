library(shiny)
#library(d3heatmap)
# library(plotly)
library(heatmaply)

shinyUI(fluidPage(
  titlePanel("CellMiner Heatmap"),
  sidebarLayout(
    sidebarPanel(
      textInput("geneList", "Gene List:", "TP53 BRAF PTEN")
    ),
    mainPanel(
      # d3heatmapOutput("heatmap")
      plotlyOutput("heatmap")
    )
  )
))

