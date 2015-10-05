library(shiny)
library(rcellminer)
library(rcellminerData)
library(jsonlite)

source(system.file("shinyTemplate", "functions.R", package="rcellminer"))
config <- jsonlite::fromJSON(system.file("shinyTemplate", "config.json", package="rcellminer"))

configDataSource <- names(config)
a <- setNames(as.character(config[[configDataSource]]$types$MolData), names(config[[configDataSource]]$types$MolData))
b <- setNames(as.character(config[[configDataSource]]$types$DrugData), names(config[[configDataSource]]$types$DrugData))
prefixChoices <- c(a, b)

#--------------------------------------------------------------------------------------------------
# Molecular data and drug activity database set up.
#--------------------------------------------------------------------------------------------------

srcContent <- list()

#----[nci60]-----------------------------------------------------------------------------
nci60Env <- new.env()
data("molData", package="rcellminerData", verbose=TRUE, envir=nci60Env)
data("drugData", package="rcellminerData", verbose=TRUE, envir=nci60Env)

nci60 <- list()

nci60MolDataTypes <- as.character(config[[configDataSource]]$types$MolData)
nci60$molPharmData <- getMolDataMatrices(getAllFeatureData(nci60Env$molData))[nci60MolDataTypes]

nci60$molPharmData[["act"]] <- exprs(getAct(nci60Env$drugData))
rownames(nci60$molPharmData$act) <- paste0("act", rownames(nci60$molPharmData$act))

nci60$drugInfo <- as(featureData(getAct(nci60Env$drugData)), "data.frame")[, c("NSC", "NAME", "MOA")]
colnames(nci60$drugInfo) <- c("ID", "NAME", "MOA")
nci60$drugInfo$ID <- as.character(nci60$drugInfo$ID)
# Handle uncommon characters
nci60$drugInfo$NAME <- iconv(enc2utf8(nci60$drugInfo$NAME), sub="byte")

stopifnot(identical(unname(removeMolDataType(rownames(nci60$molPharmData$act))), 
										nci60$drugInfo$ID))
rownames(nci60$drugInfo) <- rownames(nci60$molPharmData$act)

nci60$sampleData <- getSampleData(nci60Env$molData)
rownames(nci60$sampleData) <- nci60$sampleData$Name

nci60$tissueToSamplesMap <- getTissueToSamplesMap(getSampleData(nci60Env$molData))

nci60ColorTab <- loadNciColorSet(returnDf=TRUE)
nci60ColorTab$OncoTree1 <- nci60$sampleData$OncoTree1
nci60$tissueColorMap <- c(by(nci60ColorTab, nci60ColorTab$OncoTree1, FUN = function(x) unique(x$colors)))

srcContent[["nci60"]] <- nci60
#----------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
shinyServer(function(input, output, session) {
  
  #----[Render Application Title]------------------------------------------------------
  output$public <- renderText("CellMiner")
  #------------------------------------------------------------------------------------

  #**********************************************************************************************
	output$tabsetPanel = renderUI({		
			tabsetPanel(type="tabs", tabPanel("Plot Data", plotOutput("staticPlot")))				
	})
	#**********************************************************************************************

  # TO DO: Clean up code duplication below (xPrefixUi, yPrefixUi, etc.)
  output$xPrefixUi <- renderUI({
    switch(input$xDataset, "nci60" = selectInput("xPrefix", "x-Axis Type", choices=prefixChoices, selected="exp"))
  })

  output$yPrefixUi <- renderUI({
    switch(input$yDataset, "nci60" = selectInput("yPrefix", "y-Axis Type", choices=prefixChoices, selected="act"))
  })
  
  output$showColorTissuesUi <- renderUI({
      tissueTypes <- names(srcContent[["nci60"]][["tissueToSamplesMap"]])

      switch(input$xDataset,
             "nci60" = selectInput("showColorTissues", "Color Specific Tissues?", 
                                   choices=c("all", unique(tissueTypes)), multiple=TRUE, selected="all")
      )
  })
  
  # Static Plot
  output$staticPlot <- renderPlot({
  	if (!require(rcellminerUtils)){
  		# NOTE: There is a validate in jsonlite so the call must be distinguished
  		shiny::validate(need(input$xDataset == input$yDataset,
  									"ERROR: x and y axis data sets must be the same."))
  	}
  	shiny::validate(
  		need(validateEntry(input$xPrefix, input$xId, input$xDataset, srcContent), 
  				 paste("ERROR:", paste0(input$xPrefix, input$xId), "not found.")),
  		need(validateEntry(input$yPrefix, input$yId, input$yDataset, srcContent), 
  				 paste("ERROR:", paste0(input$yPrefix, input$yId), "not found."))
  	)

  	xData <- getFeatureData(input$xPrefix, input$xId, input$xDataset, srcContent)
  	yData <- getFeatureData(input$yPrefix, input$yId, input$yDataset, srcContent)
  	
  	makePlotStatic(xData, yData, input$showColor, input$showColorTissues, input$xDataset, srcContent)
  })
})
#-----[end of shinyServer()]-----------------------------------------------------------------------
