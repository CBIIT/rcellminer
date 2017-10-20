library(shiny)
library(d3heatmap)
library(rcellminer)
library(rcellminerElasticNet)
library(geneSetPathwayAnalysis)
library(jsonlite)
library(stringr)
library(glmnet)
library(ggplot2)
library(plotly)
#library(svglite)
#library(clusterProfiler)

#library(tooltipsterR)

if (!require(rcellminerUtils)){
	warning("rcellminerUtils package must be installed for full cross-database functionality.")
}


#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
#--------------------------------------------------------------------------------------------------
config <- jsonlite::fromJSON("config.json")
appConfig <- jsonlite::fromJSON("appConfig.json")
metaConfig <- jsonlite::fromJSON("configMeta.json")
source("modal.R")
source("appUtils.R")
source("dataLoadingFunctions.R")

#if (!is.null(appConfig$appName)){
#	appTitle <- appConfig$appName
#} else{
#	appTitle <- "CellMiner"
#}

# Construct named character vector mapping displayed data source names to
# internally used source identifiers.
dataSourceChoices <- setNames(names(config),
															vapply(config, function(x) { x[["displayName"]] }, 
																		 character(1)))

metaChoices <- setNames(names(metaConfig),
												vapply(metaConfig, function(x) { x[["displayName"]] }, 
															 character(1)))
if(!file.exists("srcContent.rds")) {
	
 for (configSrcId in names(config)){
	srcName <- config[[configSrcId]][["displayName"]]
	srcPackages <- names(config[[configSrcId]][["packages"]])
	for (pkgName in srcPackages){
		 if (!require(pkgName, character.only = TRUE)){
			  dataSourceChoices[srcName] <- NA
		  	break
	  	}
	 }
  }

 if (any(is.na(dataSourceChoices))){
	stop("Check configuration file: one or more required data source packages must be installed.")
 } 
	
} else {
	srcContent <- readRDS("srcContent.rds")
}
#--------------------------------------------------------------------------------------------------

#if("rCharts" %in% installed.packages()) {
#	options(RCHART_LIB='highcharts')	
#	library(rCharts)
#	hasRCharts <- TRUE
#} else {
#	hasRCharts <- FALSE
#}

colorSet <- loadNciColorSet(returnDf=TRUE)

###--------



#--------------------------------------------------------------------------------------------------




shinyServer(function(input, output, session) {
	#----[Reactive Variables]---------------------------------------------------------------
	# Record current input validity status, data type prefix values.
	globalReactiveValues <- reactiveValues(xPrefix = NULL, yPrefix = NULL)
	
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
	
	isPackageLoadingComplete <- reactive({
		srcContentReactive()
		
		return(TRUE)
	})

	# Provides the set of all possible feature identifiers available for plotting along
	# the x-axis of the 2-D scatter plot (and for use in pattern comparisons).
	# xIdChoices <- reactive({
	# 	srcContent <- srcContentReactive()
	# 
	# 	t1 <- lapply(srcContent[[input$xDataset]][["molPharmData"]], function(x) {
	# 		return(unname(removeMolDataType(rownames(x))))
	# 	})
	# 	l1 <- unique(unname(unlist(t1)))
	# 
	# 	t2 <- srcContent[[input$xDataset]][["drugInfo"]]
	# 	l2 <- unname(removeMolDataType(rownames(t2)))
	# 
	# 	l3 <- c(l1, l2)
	# 
	# 	return(l3)
	# })

	# Provides the set of all possible feature identifiers available for plotting along
	# the y-axis of the 2-D scatter plot.
	# yIdChoices <- reactive({
	# 	srcContent <- srcContentReactive()
	# 
	# 	t1 <- lapply(srcContent[[input$yDataset]][["molPharmData"]], function(x) {
	# 		return(unname(removeMolDataType(rownames(x))))
	# 	})
	# 	l1 <- unique(unname(unlist(t1)))
	# 
	# 	t2 <- srcContent[[input$yDataset]][["drugInfo"]]
	# 	l2 <- unname(removeMolDataType(rownames(t2)))
	# 
	# 	l3 <- c(l1, l2)
	# 
	# 	return(l3)
	# })
	
	
	# Returns all valid OncoTree types with respect to
	# the xDataset cell line OncoTree types AND user-selected tissue
	# types for inclusion or exclusion. These types need not be
	# mutually exclusive, given the nested OncoTree structure.
	# NOTE: with added data packages, always verify that matched cell
	# line OncoTree tissue type annotations are consistent across packages.
	analysisTissueTypes <- reactive({
		# Note: We want this code to re-run whenever either 
		# input$tissueSelectionMode OR  input$selectedTissues change.
		# BUT, *reactivity can be based on the selectedTissues alone*,
		# because switching the tissueSelectionMode will always trigger
		# a change to a distinct (default) selectedTissues value (e.g.,
		# "all" in the case of "Include", or "none" in the case of 
		# "Exclude"). 
		# Without the isolate() around input$tisssueSelectionMode (below),
		# this code actually runs twice upon a change to the
		# tissueSelectionMode, with the first run having an invalid
		# selectedTissues value (because the reativity relative
		# to the tissueSelectionMode appears to trigger the code below
		# before the selectedValues changes, for the first run). Then
		# the update of the selectedValues to the default value triggers
		# a second run. This behavior causes a bug-like re-drawing of
		# the 2D plot, etc.)
		tissueSelectionMode <- isolate(input$tissueSelectionMode)
		selectedTissues <- input$selectedTissues
		
		# cat("--- Entering analysisTissueTypes()", sep = "\n")
		# cat(paste0("Selection Mode: ", tissueSelectionMode), sep = "\n")
		# cat(paste0("Selected Tissues: ", selectedTissues), sep = "\n")
		
		srcContent <- srcContentReactive()
		tissueToSamplesMap <- srcContent[[input$xDataset]][["tissueToSamplesMap"]]
		tissueTypes <- names(tissueToSamplesMap)
		
		if (tissueSelectionMode == "Include"){
			if (!("all" %in% selectedTissues)){
				selectedLines <- unique(c(tissueToSamplesMap[selectedTissues], recursive = TRUE))
				# For which tissue types are ALL lines in selectedLines?
				allInSelectedLines <- vapply(tissueToSamplesMap, function(x){
					all(x %in% selectedLines)
				}, logical(1))
				tissueTypes <- names(tissueToSamplesMap[allInSelectedLines])
			}
		} else{ # tissueSelectionMode == "Exclude"
			if (!("none" %in% selectedTissues)){
				selectedLines <- unique(c(tissueToSamplesMap[selectedTissues], recursive = TRUE))
				# For which tissue types are NO lines in selectedLines?
				notInSelectedLines <- vapply(tissueToSamplesMap, function(x){
					length(intersect(x, selectedLines)) == 0
				}, logical(1))
				tissueTypes <- names(tissueToSamplesMap[notInSelectedLines])
			}
		}
		
		#cat("--- LEAVING analysisTissueTypes()", sep = "\n")
		
		return(sort(unique(tissueTypes)))
	})
	
	# Provides a data frame with columns indicating the matched cell lines between
	# the input$xDataset (column 1) and the input$yDataset (column 2).
	# matchedCellLinesTab <- reactive({
	# 	shiny::validate(need(require(rcellminerUtils),
	# 											 "ERROR: x and y axis data sets must be the same."))
	# 	matchedCellLinesTab <- getMatchedCellLines(c(input$xDataset, input$yDataset))
	# 	shiny::validate(need(nrow(matchedCellLinesTab) > 0, 
	# 											 "There are no shared cell lines between the selected datasets."))
	# 	colnames(matchedCellLinesTab) <- c("xDataset", "yDataset")
	# 	return(matchedCellLinesTab)
	# })
	
	# Provides a data frame with columns indicating the matched cell lines between
	# the input$xDataset (column 1) and the input$yDataset (column 2).
	# The cell line pairing will be updated to reflect restrictions based on:
	# --- matched cell lines across databases (if input$xDataset != input$yDataset)
	# --- user tissue type selections.
	matchedCellLinesTab <- reactive({
		srcContent <- srcContentReactive()
		analysisTissueTypes <- analysisTissueTypes()
		
		if (input$xDataset == input$yDataset){
			matchedCellLinesTab <- data.frame(
				xDataset = srcContent[[input$xDataset]]$sampleData[, "Name"],
				stringsAsFactors = FALSE
			)
			matchedCellLinesTab$yDataset <- matchedCellLinesTab$xDataset
		} else{
			shiny::validate(need(require(rcellminerUtils),
													 "ERROR: x and y axis data sets must be the same."))
			matchedCellLinesTab <- getMatchedCellLines(c(input$xDataset, input$yDataset))
			shiny::validate(need(nrow(matchedCellLinesTab) > 0, 
													 "There are no shared cell lines between the selected datasets."))
			colnames(matchedCellLinesTab) <- c("xDataset", "yDataset")
		}
		stopifnot(all(!duplicated(matchedCellLinesTab$xDataset)))
		rownames(matchedCellLinesTab) <- matchedCellLinesTab$xDataset
		
		tissueMatchedLines <- getTissueTypeSamples(analysisTissueTypes, input$xDataset, srcContent)
		tissueMatchedLines <- intersect(tissueMatchedLines, rownames(matchedCellLinesTab))
		
		shiny::validate(need(length(tissueMatchedLines) > 0, 
												 "There are no cell lines of the selected tissue type(s)."))
		
		matchedCellLinesTab <- matchedCellLinesTab[tissueMatchedLines, ]
		
		return(matchedCellLinesTab)
	})
	
	# Explanation of xData, yData reactive variables -------------------------------------------------
	# The xData and yData reactive variables provide list objects (accessed via xData() and yData()) 
	# that store the essential information about a data source feature that the application code
	# requires (e.g., to make plots, for pattern comparisons, etc.).
	# For example, if the x-axis feature is SLFN11 NCI-60 mRNA expression, the xData() 
	# list object would contain:
	# xData()$dataSource: "nci60"
	# xData()$name: "expSLFN11" (the feature identifier prepended with a data type specifier)
	# xData()$uniqName: "expSLFN11_nci60" (the above identifier appended with the data source)
	# xData()$plotLabel: "SLFN11 (exp, nci60)"
	# xData()$data: A vector of numeric feature data.
	#
	# Note: the reactive variable construction code ensures that:
	# (1) User-requested feature data is available for the specified data source,
	# (2) if the x-axis and y-axis features are derived from different data sources,
	#     the xData()$data and yData()$data vectors *have feature data from matched cell lines*.
	#
	# Through the reactivity, these variables are always 'current', reflecting the latest
	# user selections, with regard to feature identifiers, data source(s), etc.
	# ------------------------------------------------------------------------------------------------
	
	# Provides an a list object with x-axis feature-related data, including numeric data,
	# data type prefix, data source, and plot label.
	# Note that use of reactive matchedLinesTab() ensures that xData() and yData() are
	# current and coupled, reflecting matched cell lines (even when different x any y data 
	# sources are selected) of whatever tissue types are selected.
	xData <- reactive({
		shiny::validate(need(length(input$selectedTissues) > 0, "Please select tissue types."))

		xPrefix <- input$xPrefix
		if (!is.character(xPrefix)){
			xPrefix <- srcContentReactive()[[input$xDataset]][["defaultFeatureX"]]
		}
		
		xId <- getMatchedIds(xPrefix, trimws(input$xId), input$xDataset, srcContent = srcContentReactive())
		
		if (length(xId) == 0){
			shiny::validate(need(FALSE, paste("ERROR:", paste0(xPrefix, input$xId), "not found. Please use the Search IDs tab to find available IDs for each dataset.")))
		} else{
			globalReactiveValues$xPrefix <- xPrefix
			if (length(xId) > 1){
				warningMsg <- paste0("Other identifiers matching x-axis ID: ",
														 paste0(xId[-1], collapse = ", "), ".")
				showNotification(warningMsg, duration = 10, type = "message")
				xId <- xId[1]
			}
			xData <- getFeatureData(xPrefix, xId, input$xDataset, srcContent = srcContentReactive())
			
			matchedLinesTab <- matchedCellLinesTab()
			xData$data <- xData$data[matchedLinesTab[, "xDataset"]]
		}
		
		return(xData)
	})
	
	# Provides an a list object with y-axis feature-related data, including numeric data,
	# data type prefix, data source, and plot label.
	# Note that use of reactive matchedLinesTab() ensures that xData() and yData() are
	# current and coupled, reflecting matched cell lines (even when different x any y data 
	# sources are selected) of whatever tissue types are selected.
	yData <- reactive({
		shiny::validate(need(length(input$selectedTissues) > 0, "Please select tissue types."))
		
		yPrefix <- input$yPrefix
		if (!is.character(yPrefix)){
			yPrefix <- srcContentReactive()[[input$yDataset]][["defaultFeatureY"]]
		}
		
		yId <- getMatchedIds(yPrefix, trimws(input$yId), input$yDataset, srcContent = srcContentReactive())
		
		if (length(yId) == 0){
			shiny::validate(need(FALSE, paste("ERROR:", paste0(yPrefix, input$yId), "not found. Please use the Search IDs tab to find available IDs for each dataset.")))
		} else{
			globalReactiveValues$yPrefix <- yPrefix
			if (length(yId) > 1){
				warningMsg <- paste0("Other identifiers matching y-axis ID: ",
														 paste0(yId[-1], collapse = ", "), ".")
				showNotification(warningMsg, duration = 10, type = "message")
				yId <- yId[1]
			}
			yData <- getFeatureData(yPrefix, yId, input$yDataset, srcContent = srcContentReactive())
			
			matchedLinesTab <- matchedCellLinesTab()
			yData$data <- yData$data[matchedLinesTab[, "yDataset"]]
		}
		
		return(yData)
	})
	
	#----[outputs]--------------------------------------------------------------------------

  #----[Render Application Title]------------------------------------------------------
  output$public <- renderText("CellMiner")
  #------------------------------------------------------------------------------------

  #----[Render 2D Plot in 'Plot Data' Tab]---------------------------------------------
  if(require(rCharts)) {
		output$rCharts <- renderChart({
			h1 <- makePlot(xData = xData(), yData = yData(), showColor = input$showColor,
										 showColorTissues = input$showColorTissues, dataSource = input$xDataset,
										 srcContent = srcContentReactive(), dom="rCharts")
		})
  }

	# Alternative plotting
	output$rChartsAlternative <- renderPlotly({
		#-----[range check]----------------------------------------------------------
		# Note: Until a better solution can be found, these checks are needed.
		# The issue is that upon a data source change, there appears to be a moment 
		# when the xData() or yData() are updated, but the input$xAxisRange
		# or input$yAxisRange (from the sliderInput UI element) are not yet updated.
		# As such, the data value range can be out of synch with the invalidated
		# axis ranges. In the most extreme case, there are no points in the
		# specified range. The reactive code is quickly re-run with the proper
		# inputs, correcting the plot, but the error flashes briefly in a buggy 
		# looking way. 
		# Below we do a range check and quietly exit if something is amiss (knowing
		# that the reactivity will ensure that the code is re-run with a proper
		# set of inputs once thing settle down).
		#****************************************************************************
		xData <- xData()
		yData <- yData()
		shiny::validate(need(xData$uniqName != yData$uniqName, 
			"Please select distinct x and y axis variables."))
		
		xValRange <- range(xData$data, na.rm = TRUE)
		xLimits <- input$xAxisRange
		
		yValRange <- range(yData$data, na.rm = TRUE)
		yLimits <- input$yAxisRange
		
		# req(FALSE) causes immediate but quiet exit.
		req(!any(is.null(xValRange)), !any(is.null(xLimits)))
		req(!any(is.null(yValRange)), !any(is.null(yLimits)))
		req(!any(is.na(xValRange)), !any(is.na(xLimits)))
		req(!any(is.na(yValRange)), !any(is.na(yLimits)))
		
		req((xLimits[1] <= xValRange[1]) && (xValRange[2] <= xLimits[2]))
		req((yLimits[1] <= yValRange[1]) && (yValRange[2] <= yLimits[2]))
		
		cat("xAxis Limits: ", paste0(xLimits, collapse = " "), sep = "\n")
		cat("X_VAL_RANGE: ",  paste0(xValRange, collapse = " "), sep = "\n")
		
		cat("yAxis Limits: ", paste0(yLimits, collapse = " "), sep = "\n")
		cat("Y_VAL_RANGE: ",  paste0(yValRange, collapse = " "), sep = "\n")
		cat("-------------------------------------------", sep = "\n")
		#----------------------------------------------------------------------------
		
		p1 <- makePlotStatic(xData = xData, yData = yData, showColor = input$showColor, 
												 showColorTissues = input$showColorTissues, dataSource = input$xDataset, 
												 xLimVals = xLimits, yLimVals = yLimits,
												 srcContent = srcContentReactive())
		g1 <- ggplotly(p1, width=plotWidth, height=plotHeight, tooltip=tooltipCol)
		g1 <- layout(g1, margin=list(t = 75))
		g2 <- config(p = g1, collaborate=FALSE, cloud=FALSE, displaylogo=FALSE,
								 modeBarButtonsToRemove=c("select2d", "sendDataToCloud", "pan2d", "resetScale2d",
								 												 "hoverClosestCartesian", "hoverCompareCartesian",
								 												 "lasso2d", "zoomIn2d", "zoomOut2d"))
		g2
	})
	#--------------------------------------------------------------------------------------

  #----[Render Data Table in 'Download Data' Tab]----------------------------------------
  # Generate an HTML table view of the data
  output$table <- DT::renderDataTable({
  	# Column selection below is to restrict to cell line, x, y features,
  	# and tissue type information (source-provided + OncoTree).
  	dlDataTab <- getPlotData(xData = xData(), yData = yData(), showColor = input$showColor, 
  		showColorTissues = input$showColorTissues, dataSource = input$xDataset, 
  		srcContent = srcContentReactive())
  	dlDataTabCols <- c(colnames(dlDataTab)[1:4], paste0("OncoTree", 1:4))
  	if ("EMT" %in% colnames(dlDataTab)) {
  		dlDataTabCols <- c(dlDataTabCols, "EMT")
  	}
  	dlDataTab <- dlDataTab[, dlDataTabCols]
  	dlDataTab[, 2] <- round(dlDataTab[, 2], 3)
  	dlDataTab[, 3] <- round(dlDataTab[, 3], 3)

  	DT::datatable(dlDataTab, rownames=FALSE, colnames=colnames(dlDataTab),
  								filter='top', style='bootstrap', selection="none",
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
    							filter='top', style='bootstrap', selection = "none",
    							options=list(pageLength = 10))
  })
	#--------------------------------------------------------------------------------------

	#----[Render Data Table in 'Compare Patterns' Tab]-------------------------------------
	output$patternComparison <- DT::renderDataTable({
		srcContent <- srcContentReactive()
		
		if (input$patternComparisonSeed == "xPattern"){
			dat <- xData()
			pcDataset <- input$xDataset
		} else{
			dat <- yData()
			pcDataset <- input$yDataset
		}
		selectedLines <- names(dat$data)

	  if(input$patternComparisonType == "drug") {
	    results <- patternComparison(dat$data,
	    														 srcContent[[pcDataset]][["molPharmData"]][["act"]][, selectedLines])
	    results$ids <- rownames(results)
	    results$NAME <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "NAME"]

	    if ("MOA" %in% colnames(srcContent[[pcDataset]][["drugInfo"]])){
	    	results$MOA <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "MOA"]
	    	results <- results[, c("ids", "NAME", "MOA", "COR", "PVAL")]
	    	colnames(results) <- c("ID", "Name", "MOA", "Correlation", "P-Value")
	    } else{
	    	results <- results[, c("ids", "NAME", "COR", "PVAL")]
	    	colnames(results) <- c("ID", "Name", "Correlation", "P-Value")
	    }
	  } else {
	    molPharmData <- srcContent[[pcDataset]][["molPharmData"]]
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
	    
	    if (require(geneSetPathwayAnalysis)){
	    	results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[results$Gene, "SHORT_ANNOT"]
				results$Annotation[is.na(results$Annotation)] <- ""
	    }
	    
	    results$ID <- NULL
	  }
	  
	  results[, "Correlation"] <- round(results[, "Correlation"], 3)
	  results[, "P-Value"] <- signif(results[, "P-Value"], 3)
		
	  DT::datatable(results, rownames=FALSE, colnames=colnames(results),
	  							filter='top', style='bootstrap', selection = "none",
	  							options=list(lengthMenu = c(10, 25, 50, 100), pageLength = 10))
	})

	#----[Render Data Table in 'Metadata' Tab]-------------------------------------------
	output$cellLineTable <- DT::renderDataTable({
		
		configSelect <- metaConfig[[input$mdataSource]][["packages"]][[1]][["MetaData"]]
		jsonFrame <- as.data.frame(configSelect)
		
		colnames(jsonFrame) <- c("Data Type", "Description", "Units", 
														 "Platform/Assay", "PubMed Ref. ID")
		
		DT::datatable(jsonFrame, rownames=FALSE, colnames=colnames(jsonFrame),
									filter='top', style='bootstrap', selection = "none",
									options=list(pageLength = 10))
	})
	#---------------------------------------------------------------------------------------
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
										 fluidRow(
                     	column(3, selectInput("patternComparisonType", "Pattern Comparison",
                                 						choices=c("Molecular Data"="molData", "Drug Data"="drug"), 
                     												selected="molData")),
										 	column(3, selectInput("patternComparisonSeed", "With Respect to",
										 												choices=c("x-Axis Entry"="xPattern", 
										 																	"y-Axis Entry"="yPattern"), 
										 												selected="xPattern"))
										 ),
                     DT::dataTableOutput("patternComparison"))

		#if(input$hasRCharts == "TRUE") {
		if (FALSE) {
			tsPanel <- tabsetPanel(type="tabs",
									#tabPanel("Plot Data", htmlOutput("genUrl"), showOutput("rCharts", "highcharts")),
									tabPanel("Plot Data", showOutput("rCharts", "highcharts")),
									tab1, tab2, tab3
			)
		} else {
			plotPanel <- tabPanel("Plot Data", plotlyOutput("rChartsAlternative", width = plotWidth, height = plotHeight),
														br(), br(), p("Plot point tooltips provide additional information."))
			tsPanel <- tabsetPanel(plotPanel, tab1, tab2, tab3)
		}

		return(tsPanel)
	})
	
	#**********************************************************************************************
	output$metadataPanel = renderUI({
		#verbatimTextOutput("log") can be used for debugging
		#tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),
		
		#mtab1 <- tabPanel("Features",
											#DT::dataTableOutput("featTable"))
		mtab2 <- tabPanel("Cell Line Information",
											DT::dataTableOutput("cellLineTable"))
		#mtab3 <- tabPanel("Drug Information",
											#includeMarkdown("www/files/help.md"),
											#DT::dataTableOutput("drugTable"))
		
		dataFullname <- metaConfig[[input$mdataSource]][["fullName"]]
		
		tabsetPanel(type="pills",
								tabPanel(dataFullname, tags$hr(),
													mtab2
								)
		)
	})
	#**************************************************************************************
	output$sourceLink <- renderUI({
		
		urlString <- metaConfig[[input$mdataSource]][["url"]]
		sourceName <- metaConfig[[input$mdataSource]][["displayName"]]
		visibleText <- paste("Select here to learn more about ", sourceName, sep="")
	
		a(visibleText, href=paste(urlString), target = "_blank")
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
    	df <- getPlotData(xData = xData(), yData = yData(), showColor = input$showColor, 
    		showColorTissues = input$showColorTissues, dataSource = input$xDataset, 
    		srcContent = srcContentReactive())

    	# Column selection below is to restrict to cell line, x, y features,
    	# and tissue type information (source-provided + OncoTree).
    	dfCols <- c(colnames(df)[1:4], paste0("OncoTree", 1:4))
    	if ("EMT" %in% colnames(df)) {
    		dfCols <- c(dfCols, "EMT")
    	}
    	df <- df[, dfCols]

    	write.table(df, file, quote=FALSE, row.names=FALSE, sep="\t")
    }
  )

  output$xPrefixUi <- renderUI({
  	srcContent <- srcContentReactive()
  	
  	# The last selected (data type) prefix is recorded in 
  	# globalReactiveValues$xPrefix whenever xData() is updated. When the data set 
  	# is changed, we try to use this same data type prefix, if it is available.
  	prefixChoices <- srcContent[[input$xDataset]][["featurePrefixes"]]
  	selectedPrefix <- globalReactiveValues$xPrefix
  	if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
  		selectedPrefix <- srcContent[[input$xDataset]][["defaultFeatureX"]]
  	}
  	selectInput("xPrefix", "x-Axis Type", choices = prefixChoices, selected = selectedPrefix)
  })

  output$yPrefixUi <- renderUI({
  	srcContent <- srcContentReactive()
  	prefixChoices <- srcContent[[input$yDataset]][["featurePrefixes"]]
  	selectedPrefix <- globalReactiveValues$yPrefix
  	if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
  		selectedPrefix <- srcContent[[input$yDataset]][["defaultFeatureY"]]
  	}
  	selectInput("yPrefix", "y-Axis Type", choices = prefixChoices, selected = selectedPrefix)
  })
  
  output$xAxisRangeUi <- renderUI({
  	srcContent <- srcContentReactive()
  	
  	# Note: req() ensures values are available or 'truthy' (not NULL, "", FALSE, empty, etc.),
  	# returning the value if so; otherwise the operation is stopped with a silent exception.
  	# The idea is to exit quietly if inputs are momentarily in an invalid state, as might
  	# occur when the app is first loading, etc.
  	valRange <- srcContent[[req(input$xDataset)]][["featureValRanges"]][[req(input$xPrefix)]]
  	
  	xData <- NULL
  	try(xData <- xData())
  	if (is.null(xData)){
  		xInitSliderVals <- valRange
  	} else{
  		xDataRange <- range(xData$data, na.rm = TRUE)
  		delta <- max((0.05 * (xDataRange[2] - xDataRange[1])), 0.1)
  		xInitSliderVals <- c((xDataRange[1] - delta), (xDataRange[2] + delta))
  	}
  	
  	sliderInput("xAxisRange", "x-Axis Range", 
  							min = valRange[1], max = valRange[2], value = xInitSliderVals, step = 0.5)
  })
  
  output$yAxisRangeUi <- renderUI({
  	srcContent <- srcContentReactive()
  	
 		# Note: see comment in output#xAxisRangeUi explaining the use of req().
  	valRange <- srcContent[[req(input$yDataset)]][["featureValRanges"]][[req(input$yPrefix)]]
  	
  	yData <- NULL
  	try(yData <- yData())
  	if (is.null(yData)){
  		yInitSliderVals <- valRange
  	} else{
  		yDataRange <- range(yData$data, na.rm = TRUE)
  		delta <- max((0.05 * (yDataRange[2] - yDataRange[1])), 0.1)
  		yInitSliderVals <- c((yDataRange[1] - delta), (yDataRange[2] + delta))
  	}
  	
  	sliderInput("yAxisRange", "y-Axis Range", 
  							min = valRange[1], max = valRange[2], value = yInitSliderVals, step = 0.5)
  })

  # output$xIdUi <- renderUI({
  # 	updateSelectizeInput(session, inputId='xId', choices=xIdChoices(), selected="SLFN11", server=TRUE)
  # 	selectizeInput('xId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5))
  # })
  # 
  # output$yIdUi <- renderUI({
  # 	updateSelectizeInput(session, inputId='yId', choices=yIdChoices(), selected="94600", server=TRUE)
  # 	selectizeInput('yId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5))
  # })
  
  output$selectTissuesUi <- renderUI({
  	srcContent <- srcContentReactive()
  	tissueToSamplesMap <- srcContent[[input$xDataset]][["tissueToSamplesMap"]]
  	tissueTypes <- sort(unique(names(tissueToSamplesMap)))
  	
  	if (input$tissueSelectionMode == "Include"){
  		selectInput("selectedTissues", label = NULL, choices=c("all", tissueTypes),
  								multiple=TRUE, selected="all")
  	} else{ # input$tissueSelectionMode == "Exclude"
  		selectInput("selectedTissues", label = NULL, choices=c("none", tissueTypes),
  								multiple=TRUE, selected="none")
  	}
  })

  output$showColorTissuesUi <- renderUI({
  	#tissueChoices <- analysisTissueTypes()
  	tissueChoices <- getSampleSetTissueTypes(
  		sampleSet = rownames(req(matchedCellLinesTab())), 
  		dataSource = input$xDataset, 
  		srcContent = srcContentReactive()
  		) 
  	selectInput("showColorTissues", "Tissues to Color",
  							choices = tissueChoices, multiple = TRUE)
	})
  
  output$ipAddress <- renderText({
  	# debug
  	text <- readLines("http://api.ipify.org")
  })

  #----[observers]-----------------------------------------------------------------------

  # Observe reactive variable and send message to Javascript code
  observe({
  	if(isPackageLoadingComplete()) {
  		session$sendCustomMessage(type='showLoading', list(show=FALSE))
  	}
  })
	
  #-----[NavBar Tab Server Code]---------------------------------------------------------
  rm <- callModule(regressionModels, "rm", srcContentReactive = srcContentReactive,
  								 appConfig = appConfig)

})
#-----[end of shinyServer()]-----------------------------------------------------------------------
