library(shiny)
library(rcellminer)
library(jsonlite)

# Note: The jsonlite package has a validate() function whose name clashes with a
# validate() function provided by the shiny package.
# shiny::validate() must be used in the following code.

#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND CONSTRUCT APP (MOLECULAR + DRUG) DATA.
#--------------------------------------------------------------------------------------------------
#config <- jsonlite::fromJSON(system.file("shinyComparePlots", "config.json", package="rcellminer"))
config <- jsonlite::fromJSON("config.json")
source("appUtils.R")
source("dataLoadingFunctions.R")
#--------------------------------------------------------------------------------------------------

if(file.exists("srcContent.rds")) {
        srcContent <- readRDS("srcContent.rds")
}

shinyServer(function(input, output, session) {
	#----[Reactive Variables]---------------------------------------------------------------
	isPackageLoadingComplete <- reactive({
		srcContentReactive()

		return(TRUE)
	})

	# Provides a srcContent (list-based) data structure containing all molecular profiling
	# drug response, and feature/sample annotation data required for the application 
	# (for data sources specified in the config.json file).
	srcContentReactive <- reactive({
		if(!exists("srcContent")) {
			srcContent <- lapply(config, loadSourceContent)
			isLoadedSrc <- vapply(srcContent, function(x) { !is.null(x) }, logical(1))
			if (any(!isLoadedSrc)){
				srcContent <- srcContent[isLoadedSrc]
			}

			# For NCI-60, replace default color map to use CellMiner tissue type colors.
			nci60ColorTab <- loadNciColorSet(returnDf=TRUE)
			nci60ColorTab$OncoTree1 <- srcContent$nci60$sampleData$OncoTree1
			srcContent$nci60$tissueColorMap <- c(by(nci60ColorTab, nci60ColorTab$OncoTree1,
																							FUN = function(x) unique(x$colors)))
		}
		
		return(srcContent)
	})

	# Provides the set of all possible feature identifiers available for plotting along
	# the x-axis of the 2-D scatter plot (and for use in pattern comparisons).
	xIdChoices <- reactive({
		srcContent <- srcContentReactive()

		t1 <- lapply(srcContent[[input$xDataset]][["molPharmData"]], function(x) {
			return(unname(removeMolDataType(rownames(x))))
		})
		l1 <- unique(unname(unlist(t1)))

		t2 <- srcContent[[input$xDataset]][["drugInfo"]]
		l2 <- unname(removeMolDataType(rownames(t2)))

		l3 <- c(l1, l2)

		return(l3)
	})

	# Provides the set of all possible feature identifiers available for plotting along
	# the y-axis of the 2-D scatter plot.
	yIdChoices <- reactive({
		srcContent <- srcContentReactive()

		t1 <- lapply(srcContent[[input$yDataset]][["molPharmData"]], function(x) {
			return(unname(removeMolDataType(rownames(x))))
		})
		l1 <- unique(unname(unlist(t1)))

		t2 <- srcContent[[input$yDataset]][["drugInfo"]]
		l2 <- unname(removeMolDataType(rownames(t2)))

		l3 <- c(l1, l2)

		return(l3)
	})
	
	# Provides an a list object with x-axis feature-related data, including numeric data,
	# data type prefix, data source, and plot label.
	xData <- reactive({
		shiny::validate(need(validateEntry(input$xPrefix, input$xId, input$xDataset,
												 srcContent = srcContentReactive()),
												 paste("ERROR:", paste0(input$xPrefix, input$xId), "not found.")))
		xData <- getFeatureData(input$xPrefix, input$xId, input$xDataset, 
														srcContent = srcContentReactive())
		return(xData)
	})
	
	# Provides an a list object with y-axis feature-related data, including numeric data,
	# data type prefix, data source, and plot label.
	yData <- reactive({
		shiny::validate(need(validateEntry(input$yPrefix, input$yId, input$yDataset,
												 srcContent = srcContentReactive()),
												 paste("ERROR:", paste0(input$yPrefix, input$yId), "not found.")))
		yData <- getFeatureData(input$yPrefix, input$yId, input$yDataset, 
														srcContent = srcContentReactive())
		return(yData)
	})

	#----[outputs]--------------------------------------------------------------------------

  #----[Render Application Title]------------------------------------------------------
  output$public <- renderText("CellMiner")
  #------------------------------------------------------------------------------------

  #----[Render 2D Plot in 'Plot Data' Tab]---------------------------------------------
  if(require(rCharts)) {
		output$rCharts <- renderChart({
			if (input$selectedTissuesOnly){
				shiny::validate(need(length(input$showColorTissues) > 0, "Please select tissue types."))
			}
			if (!require(rcellminerUtils)){
				shiny::validate(need(input$xDataset == input$yDataset,
							   "ERROR: x and y axis data sets must be the same."))
			}

			h1 <- makePlot(xData = xData(), yData = yData(), showColor = input$showColor,
										 showColorTissues = input$showColorTissues, dataSource = input$xDataset,
										 selectedTissuesOnly = input$selectedTissuesOnly, 
										 srcContent = srcContentReactive())
		})
  }

	# Alternative plotting
	output$rChartsAlternative <- renderPlot({
		if (input$selectedTissuesOnly){
			shiny::validate(need(length(input$showColorTissues) > 0, "Please select tissue types."))
		}
		if (!require(rcellminerUtils)){
			shiny::validate(need(input$xDataset == input$yDataset,
										"ERROR: x and y axis data sets must be the same."))
		}

		makePlotStatic(xData(), yData(), input$showColor, input$showColorTissues, input$xDataset,
									 input$selectedTissuesOnly, srcContent = srcContentReactive())
	})
	#--------------------------------------------------------------------------------------

  #----[Render Data Table in 'Download Data' Tab]----------------------------------------
  # Generate an HTML table view of the data
  output$table <- DT::renderDataTable({
  	if (input$selectedTissuesOnly){
  		shiny::validate(need(length(input$showColorTissues) > 0, "Please select tissue types."))
  	}
  	if (!require(rcellminerUtils)){
  		shiny::validate(need(input$xDataset == input$yDataset, "ERROR: x and y axis data sets must be the same."))
  	}

  	# Column selection below is to restrict to cell line, x, y features,
  	# and tissue type information (source-provided + OncoTree).
  	dlDataTab <- getPlotData(xData(), yData(), input$showColor, input$showColorTissues,
  													 input$xDataset, input$selectedTissuesOnly, srcContent = srcContentReactive())
  	dlDataTabCols <- c(colnames(dlDataTab)[1:4], paste0("OncoTree", 1:4))
  	dlDataTab <- dlDataTab[, dlDataTabCols]

  	DT::datatable(dlDataTab, rownames=FALSE, colnames=colnames(dlDataTab),
  								filter='top', style='bootstrap',
  								options=list(pageLength = nrow(dlDataTab)))
  })
	#--------------------------------------------------------------------------------------

  #----[Render Data Table in 'Search IDs' Tab]-------------------------------------------
  # Generate an HTML table view of the data
  # Note: Searchable data is derived from the x-axis data source.
	output$ids <- DT::renderDataTable({
		srcContent <- srcContentReactive()
    drugIds   <- srcContent[[input$xDataset]][["drugInfo"]][, "ID"]
    drugNames <- srcContent[[input$xDataset]][["drugInfo"]][, "NAME"]
    moaNames  <- srcContent[[input$xDataset]][["drugInfo"]][, "MOA"]

    results <- data.frame(availableIds = drugIds, idNames=drugNames, properties=moaNames,
                          stringsAsFactors=FALSE)

    # Make molecular data data.frame
    molPharmData <- srcContent[[input$xDataset]][["molPharmData"]]
    molData <- molPharmData[setdiff(names(molPharmData), "act")]
    molDataIds <- as.vector(unlist(lapply(molData, function(x) { rownames(x) })))
    molDataNames <- rep("", length(molDataIds))
    moaNames <- rep("", length(molDataIds))

    tmp <- data.frame(availableIds=molDataIds, idNames=molDataNames, properties=moaNames,
                      stringsAsFactors=FALSE)

    # Join data.frames
    results <- rbind(results, tmp)

    # Reverse Order
    results <- results[rev(rownames(results)),]
    colnames(results) <- c("ID", "Name", "Drug MOA")

    DT::datatable(results, rownames=FALSE, colnames=colnames(results),
    							filter='top', style='bootstrap',
    							options=list(pageLength = 10))
  })
	#--------------------------------------------------------------------------------------

	#----[Render Data Table in 'Compare Patterns' Tab]-------------------------------------
	output$patternComparison <- DT::renderDataTable({
		srcContent <- srcContentReactive()
		if (input$selectedTissuesOnly){
			shiny::validate(need(length(input$showColorTissues) > 0, "Please select tissue types."))
		}
		if (!require(rcellminerUtils)){
			shiny::validate(need(input$xDataset == input$yDataset, "ERROR: x and y axis data sets must be the same."))
		}
	  # shiny::validate(need(validateEntry(input$xPrefix, input$xId, input$xDataset,
	  # 																	 srcContent = srcContentReactive()),
	  # 										 "ERROR: x-Axis Entry Not Found."))
	  # 
	  # dat <- getFeatureData(input$xPrefix, input$xId, input$xDataset, srcContent = srcContentReactive())
		dat <- xData()
		if (require(rcellminerUtils)){
			matchedLinesTab <- getMatchedCellLines(c(input$xDataset, input$yDataset))
			dat$data <- dat$data[matchedLinesTab[, 1]]
		}

		if (input$selectedTissuesOnly){
			if ((length(input$showColorTissues) > 0) && (!("all" %in% input$showColorTissues))){
				matchedSamples <- getTissueTypeSamples(input$showColorTissues, input$xDataset, srcContent)
				dat$data <- dat$data[intersect(matchedSamples, names(dat$data))]
			}
		}

	  selectedLines <- names(dat$data)

	  if(input$patternComparisonType == "drug") {
	    results <- patternComparison(dat$data,
	    														 srcContent[[input$xDataset]][["molPharmData"]][["act"]][, selectedLines])
	    results$ids <- rownames(results)
	    results$NAME <- srcContent[[input$xDataset]][["drugInfo"]][rownames(results), "NAME"]

	    if ("MOA" %in% colnames(srcContent[[input$xDataset]][["drugInfo"]])){
	    	results$MOA <- srcContent[[input$xDataset]][["drugInfo"]][rownames(results), "MOA"]
	    	results <- results[, c("ids", "NAME", "MOA", "COR", "PVAL")]
	    	colnames(results) <- c("ID", "Name", "MOA", "Correlation", "P-Value")
	    } else{
	    	results <- results[, c("ids", "NAME", "COR", "PVAL")]
	    	colnames(results) <- c("ID", "Name", "Correlation", "P-Value")
	    }
	  } else {
	    molPharmData <- srcContent[[input$xDataset]][["molPharmData"]]
	    molData <- molPharmData[setdiff(names(molPharmData), "act")]
	    molData <- lapply(molData, function(X) X[, selectedLines])
	    results <- patternComparison(dat$data, molData)
	    results$ids <- rownames(results)

	    results$molDataType <- getMolDataType(results$ids)
	    results$gene <- removeMolDataType(results$ids)

	    # Reorder columns
	    results <- results[, c("ids", "molDataType", "gene", "COR", "PVAL")]
	    colnames(results) <- c("ID", "Data Type", "Gene", "Correlation", "P-Value")

	    if (require(rcellminerUtils)){
	    	chromLocs <- character(nrow(results))
	    	haveLoc <- results$Gene %in% names(geneToChromBand)
	    	chromLocs[haveLoc] <- geneToChromBand[results$Gene[haveLoc]]

	    	results$Location <- chromLocs
	    	results <- results[, c("ID", "Data Type", "Gene", "Location", "Correlation", "P-Value")]
	    }
	  }

	  DT::datatable(results, rownames=FALSE, colnames=colnames(results),
	  							filter='top', style='bootstrap',
	  							options=list(lengthMenu = c(10, 25, 50, 100), pageLength = 10))
	})
	#--------------------------------------------------------------------------------------

	output$log <- renderText({
			paste(names(input), collapse=" ")
			query <- parseQueryString(session$clientData$url_search)
			sapply(names(input), function(x) { query[[x]] <- input[[x]] })
		})

# 	output$genUrl <- renderText({
# 		#query <- parseQueryString(session$clientData$url_search)
# 		query <- list()
#
# 		# Get current input values and combine with
# 		tmp <- sapply(names(input), function(x) { query[[x]] <- input[[x]] })
# 		query <- c(query, tmp)
#
# 		paramStr <- paste(sapply(names(query), function(x) { paste0(x, "=", query[[x]]) }), collapse="&")
#
# 		urlStr <- paste(session$clientData$url_protocol, "//",
# 										session$clientData$url_hostname, ":",
# 										session$clientData$url_port,
# 										session$clientData$url_pathname,
# 										"?", paramStr,
# 										sep="")
#
# 		paste0(a(paste("Shareable Link (Right-click then 'Copy' to share)"), href=urlStr), hr())
#
# 		#paste("Shareable Link:", urlStr)
# 	})

  #**********************************************************************************************
	output$tabsetPanel = renderUI({
		#verbatimTextOutput("log") can be used for debugging
		#tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),

		tab1 <- tabPanel("Download Data",
                     downloadLink("downloadData", "Download Data as Tab-Delimited File"),
                     DT::dataTableOutput("table"))
		tab2 <- tabPanel("Search IDs",
                     includeMarkdown("www/files/help.md"),
                     DT::dataTableOutput("ids"))
		tab3 <- tabPanel("Compare Patterns",
										 includeMarkdown("www/files/help.md"),
                     selectInput("patternComparisonType", "Pattern Comparison to x-Axis Entry",
                                 choices=c("Molecular Data"="molData", "Drug Data"="drug"), selected="molData"),
                     DT::dataTableOutput("patternComparison"))

		if(input$hasRCharts == "TRUE") {
			tabsetPanel(type="tabs",
									#tabPanel("Plot Data", htmlOutput("genUrl"), showOutput("rCharts", "highcharts")),
									tabPanel("Plot Data", showOutput("rCharts", "highcharts")),
									tab1, tab2, tab3
			)
		} else {
			tabsetPanel(type="tabs",
									#tabPanel("Plot Data", htmlOutput("genUrl"), plotOutput("rChartsAlternative", width=600, height=600)),
									tabPanel("Plot Data", plotOutput("rChartsAlternative", width=600, height=600)),
									tab1, tab2, tab3
			)
		}
	})
	#**********************************************************************************************

  output$downloadData <- downloadHandler(
    filename = function() {
    	query <- parseQueryString(session$clientData$url_search)

    	if("filename" %in% names(query)) {
    		filename <- query[["filename"]]
    	} else {
    		filename <- "dataset"
    	}

    	if("extension" %in% names(query)) {
    		extension <- query[["extension"]]
    	} else {
    		extension <- "txt"
    	}

    	paste(filename, extension, sep=".")
    },
    content = function(file) {
    	df <- getPlotData(xData(), yData(), input$showColor, input$showColorTissues,
    										input$xDataset, input$selectedTissuesOnly, srcContent = srcContentReactive())

    	# Column selection below is to restrict to cell line, x, y features,
    	# and tissue type information (source-provided + OncoTree).
    	dfCols <- c(colnames(df)[1:4], paste0("OncoTree", 1:4))
    	df <- df[, dfCols]

    	write.table(df, file, quote=FALSE, row.names=FALSE, sep="\t")
    }
  )

  output$xPrefixUi <- renderUI({
  	srcContent <- srcContentReactive()
  	selectInput("xPrefix", "x-Axis Type",
  							choices = srcContent[[input$xDataset]][["featurePrefixes"]],
  							selected = srcContent[[input$xDataset]][["defaultFeatureX"]])
  })

  output$yPrefixUi <- renderUI({
  	srcContent <- srcContentReactive()
  	selectInput("yPrefix", "y-Axis Type",
  							choices = srcContent[[input$yDataset]][["featurePrefixes"]],
  							selected = srcContent[[input$yDataset]][["defaultFeatureY"]])
  })

  output$xIdUi <- renderUI({
  	updateSelectizeInput(session, inputId='xId', choices=xIdChoices(), selected="SLFN11", server=TRUE)
  	selectizeInput('xId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5))
  })

  output$yIdUi <- renderUI({
  	updateSelectizeInput(session, inputId='yId', choices=yIdChoices(), selected="94600", server=TRUE)
  	selectizeInput('yId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5))
  })

  output$showColorTissuesUi <- renderUI({
  	srcContent <- srcContentReactive()
  	tissueTypes <- names(srcContent[[input$xDataset]][["tissueToSamplesMap"]])
  	selectInput("showColorTissues", "Color Specific Tissues?",
  							choices=c("all", unique(tissueTypes)), multiple=TRUE, selected="all")
	})

  #----[observers]-----------------------------------------------------------------------

  # Observe reactive variable and send message to Javascript code
  observe({
  	if(isPackageLoadingComplete()) {
  		session$sendCustomMessage(type='showLoading', list(show=FALSE))
  	}
  })

})
#-----[end of shinyServer()]-----------------------------------------------------------------------
