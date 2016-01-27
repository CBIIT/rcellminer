library(shiny)
library(rcellminer)
library(d3heatmap)

shinyServer(
  function(input, output){
    output$heatmap <- renderD3heatmap({
      genes <- unlist(strsplit(input$geneList, " "))
      expData <- getAllFeatureData(rcellminerData::molData)[["exp"]]
      d3heatmap(expData[genes, 1:20], scale="column", colors="YlOrRd")
    })
  }
)
