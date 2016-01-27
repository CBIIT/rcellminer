library(shiny)
library(d3heatmap)

shinyUI(fluidPage(
  titlePanel("CellMiner Heatmap"),
  sidebarLayout(
    sidebarPanel(
      textInput("geneList", "Gene List:", "TP53 BRAF PTEN")
    ),
    mainPanel(
      d3heatmapOutput("heatmap")
    )
  )
))

