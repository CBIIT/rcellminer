library(shiny)
library(rcellminer)
#library(d3heatmap)
library(heatmaply)
library(gplots)

shinyServer(
  function(input, output){
    ## output$heatmap <- renderD3heatmap({
   output$heatmap <- renderPlotly({	
      genes <- unlist(strsplit(input$geneList, " "))
      expData <- getAllFeatureData(rcellminerData::molData)[["exp"]]
      # d3heatmap(expData[genes, 1:20], scale="column", colors="YlOrRd")
      heatmaply(expData[genes, 1:20], scale="row", grid_gap = 0.5, cellnote_color="black", colors= colorpanel(75,low="blue",mid="white",high="red"),fontsize_col = 8 ,fontsize_row = 8)
    })
  }
)
