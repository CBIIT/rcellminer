library(shiny)
library(rcellminer)
library(rcellminerData)
library(jsonlite)

source(system.file("shinyTemplate", "functions.R", package="rcellminer"))
config <- jsonlite::fromJSON(system.file("shinyTemplate", "config.json", package="rcellminer"))

# configDataSource <- "nci60"
# a <- setNames(as.character(config[[configDataSource]]$types$MolData), names(config[[configDataSource]]$types$MolData))
# b <- setNames(as.character(config[[configDataSource]]$types$DrugData), names(config[[configDataSource]]$types$DrugData))
# prefixChoices <- c(a, b)

#--------------------------------------------------------------------------------------------------
# Molecular data and drug activity database set up.
#--------------------------------------------------------------------------------------------------

srcContent <- list()

#----[nci60]-----------------------------------------------------------------------------
nci60Env <- new.env()
data("molData", package="rcellminerData", verbose=TRUE, envir=nci60Env)
data("drugData", package="rcellminerData", verbose=TRUE, envir=nci60Env)

nci60 <- list()

nci60MolDataTypes <- c("exp", "xai", "mut", "cop", "mir", "mda", "pro")
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
			tabsetPanel(type="tabs", tabPanel("Plot Data", plotOutput("formulaCorPlot"), verbatimTextOutput("formulaCorSummary")))				
	})
	#**********************************************************************************************

  # TO DO: Clean up code duplication below (xPrefixUi, yPrefixUi, etc.)
  output$xPrefixUi <- renderUI({
  	prefixChoices <- c("Expression (Z-Score)"="exp", "Expression (Avg. log2 Int.)"="xai",
  										 "Mutations"="mut", "Copy Number"="cop", 
  										 "Drug"="act", "MicroRNA"="mir", "Metadata"="mda", 
  										 "Protein (RPLA)"="pro")
    switch(input$xDataset, "nci60" = selectInput("xPrefix", "x-Axis Type", choices=prefixChoices, selected="exp"))
  })

  output$yPrefixUi <- renderUI({
  	prefixChoices <- c("Expression (Z-Score)"="exp", "Expression (Avg. log2 Int.)"="xai",
  										 "Mutations"="mut", "Copy Number"="cop", 
  										 "Drug"="act", "MicroRNA"="mir", "Metadata"="mda", 
  										 "Protein (RPLA)"="pro")
    switch(input$yDataset, "nci60" = selectInput("yPrefix", "y-Axis Type", choices=prefixChoices, selected="act"))
  })
  
  output$showColorTissuesUi <- renderUI({
      tissueTypes <- names(srcContent[["nci60"]][["tissueToSamplesMap"]])

      switch(input$xDataset,
             "nci60" = selectInput("showColorTissues", "Color Specific Tissues?", 
                                   choices=c("all", unique(tissueTypes)), multiple=TRUE, selected="all")
      )
  })
  
  output$formulaCorPlot <- renderPlot({
    shiny::validate(need(grepl("~", input$xId), "ERROR: Incorrect formula."))

    eq <- input$xId
    
    # Construct data.frame 
    terms <- unlist(strsplit(eq, '[~+]'))
    
    tmp <- NULL 
    
    for(term in terms) {
      molDataType <- getMolDataType(term)
      tmp <- rbind(tmp, nci60$molPharmData[[molDataType]][term,])
    }
    
    rownames(tmp) <- terms
    df <- as.data.frame(t(tmp))
    
    # Do regression 
    fit <- lm(as.formula(eq), data=df)
    
    # Get prediction 
    pred <- predict(fit)
    
    # Get response variable 
    molDataType <- getMolDataType(terms[1])
    y <- nci60$molPharmData[[molDataType]][terms[1],]
    
    # Plot only matching entries 
    samples <- intersect(names(pred), names(y))
    y <- y[samples]
    
    # Add regression between predicted and observed 
    fit2 <- lm(pred ~ y)
    
    # Plot predicted versus observed
    title <- paste0("Formula: ", as.character(eq), "; R-squared: ", round(summary(fit2)$r.squared, 3))
    plot(y, pred, xlim=range(c(y,pred)), ylim=range(c(y,pred)), xlab="observed", ylab="predicted",  main=title)
    
    # Add regression line
    abline(fit2, lwd=2)
  })
  
  output$formulaCorSummary <- renderPrint({
    shiny::validate(need(grepl("~", input$xId), "ERROR: Incorrect formula."))
    
    eq <- input$xId
    
    # Construct data.frame 
    terms <- unlist(strsplit(eq, '[~+]'))
    
    tmp <- NULL 
    
    for(term in terms) {
      molDataType <- getMolDataType(term)
      tmp <- rbind(tmp, nci60$molPharmData[[molDataType]][term,])
    }
    
    rownames(tmp) <- terms
    df <- as.data.frame(t(tmp))
    
    # Do regression 
    fit <- lm(as.formula(eq), data=df)
    summary(fit)
  })
})
#-----[end of shinyServer()]-----------------------------------------------------------------------
