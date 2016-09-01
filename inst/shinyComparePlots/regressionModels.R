#-----[NavBar Tab: Regression Models (UI code)]----------------------------------------------------
regressionModelsInput <- function(id, dataSourceChoices) {
	# Create a namespace function using the provided id
	ns <- NS(id)
	tabPanel("Regression Models",
					 fluidPage(
					 	sidebarLayout(
					 		sidebarPanel(
					 			width=3, 
					 			tags$div(
					 				id="rm_input_container", 
					 				selectInput(ns("dataset"), "Dataset", choices=dataSourceChoices, selected = "nci60"),
					 				uiOutput(ns("responseDataTypeUi")),
					 				textInput(ns("responseId"), "Response ID: (Case-Sensitive, e.g., 609699)", "609699"),
					 				uiOutput(ns("predDataTypesUi")),
					 				textInput(ns("predIds"), "Predictor IDS: (Case-Sensitive, e.g. SLFN11 JAG1)", "SLFN11 JAG1"),
					 				selectInput(ns("algorithm"), "Algorithm", 
					 										choices=c("Linear Regression", "Supervised Principal Components"), 
					 										selected = "Linear Regression")
					 			)
					 		),
					 		mainPanel(
					 			uiOutput(ns('tabsetPanel'))
					 		)
					 	)
					 )						 
	)
}
#-----[NavBar Tab: Regression Models (Server code)]------------------------------------------------
regressionModels <- function(input, output, session, srcContentReactive) {

	#----[Reactive Variables]--------------------------------------------------------------
	
	# Provides an a list object with response variable-related data, including numeric data,
	# data type prefix, data source, and plot label.
	# Note: Observations with missing values (in predictor or reponse variables) are removed.
	# Note: The numeric data always matches what is provided within the inputData() data frame.
	rmResponseData <- reactive({
		dataTab <- inputData() # Validation already done.
		
		yData <- list()
		yData$name <- colnames(dataTab)[2] # 1st column: cell line name, 2nd column: response variable.
		yData$data <- setNames(dataTab[, yData$name], rownames(dataTab))
		yData$plotLabel <- yData$name
		yData$uniqName  <- yData$name
		yData$dataSource <- strsplit(yData$name, split = "_")[[1]][2]

		return(yData)
	})
	
	# Returns a data frame with data suitable for predictive modeling (updated according
	# to current user selections). The first column indicates the cell line, with
	# subsequent columns containing response and predictor variables.
	# Note: Observations with missing values (in predictor or reponse variables) are removed.
	inputData <- reactive({
		shiny::validate(need(validateEntry(input$responseDataType, 
																			 input$responseId, input$dataset,
																			 srcContent = srcContentReactive()),
			paste("ERROR:", paste0(input$responseDataType, input$responseId), "not found.")))
		yData <- getFeatureData(input$responseDataType, input$responseId, input$dataset, 
														srcContent = srcContentReactive())
		yData$data <- na.exclude(yData$data)
		
		dataTab <- data.frame(CellLine = names(yData$data), stringsAsFactors = FALSE)
		rownames(dataTab) <- dataTab$CellLine
		dataTab[, yData$uniqName] <- yData$data
		
		predIds <- stringr::str_split(stringr::str_trim(input$predIds), pattern = "\\s+")[[1]] 
		for (id in predIds){
			for (dataType in input$predDataTypes){
				if (validateEntry(dataType, id, input$dataset, srcContentReactive())){
					xData <- getFeatureData(dataType, id, input$dataset, srcContentReactive())
					xData$data <- xData$data[names(yData$data)] # Match lines w/non-NA response data.
					dataTab[, xData$uniqName] <- xData$data
				} else{
					warning(paste0(dataType, id), " not found.")
				}
			}
		}
		
		dataTab <- na.exclude(dataTab)
		
		shiny::validate(need(ncol(dataTab) > 2,
												 paste("ERROR: No data for specified predictors.")))
		
		shiny::validate(need(nrow(dataTab) > 0,
												 paste("ERROR: All cell lines have missing response or predictor data.")))
		
		return(dataTab)
	})
	
	# Returns a list with (at least) the following elements:
	# predictorWts: a named vector indicating the intercept term and predictor coefficients
	#  in the regression model.
	# algorithm: the name of the algorithm.
	# predictedResponse: a vector of predicted reponse values.
	# cvPredictedResponse: a vector of cross-validation predicted response values (can be
	#   empty if cross validation is not supported).
	# techDetails: an object that can be passed to a suitable summary() function to provide
	#   a printed summary of technical results.
	algoResults <- reactive({
		dataTab <- inputData()
		rmAlgoResults <- list()
		
		# Note: refactor this function to gather and validate parameters, which are then passed
		# to a specialized implementation function that returns a standard format algoResults
		# list object.
		if (input$algorithm == "Linear Regression"){
			lmData <- dataTab[, -1] # First column has cell line names
			# Note: will have to handle issues caused by predictor names with spaces 
			# or other characters that cannot be used within a formula.
			lmFormula <- as.formula(paste0(colnames(lmData)[1], " ~ ."))
			lmFit <- lm(lmFormula, lmData)
			lmCvFit <- rcellminerElasticNet::getLmCvFit(X = as.matrix(lmData[, -1, drop = FALSE]), 
																									y = lmData[, 1, drop = TRUE], nFolds = 10, nRepeats = 1)
			# -----[assemble results]---------------------------------------------------
			rmAlgoResults$algorithm <- "Linear Regression"
			rmAlgoResults$predictorWts <- coef(lmFit)
			rmAlgoResults$predictedResponse <- predict(lmFit)
			rmAlgoResults$cvPredictedResponse <- lmCvFit$cvPred
			rmAlgoResults$techDetails <- lmFit
			
			stopifnot(identical(names(rmAlgoResults$predictedResponse), rownames(lmData)))
			stopifnot(identical(names(rmAlgoResults$cvPredictedResponse), rownames(lmData)))
			# --------------------------------------------------------------------------
		} else{
			shiny::validate(FALSE, paste("ERROR: Algorithm not available."))
		}
		
		return(rmAlgoResults)
	})
	
	# Returns a data frame with partial correlation-based pattern comparison results.
	parCorPatternCompResults <- eventReactive(input$computeParCors, {
		shiny::validate(need(length(input$pcGeneSets) > 0,
												 "Please select one or more gene sets."))
		srcContent <- srcContentReactive()
		dataTab <- inputData()
		responseData <- rmResponseData()
		
		responseVec <- responseData$data
		currentPredictorData <- t(as.matrix(dataTab[, c(-1, -2), drop = FALSE]))
		
		comparisonData <- srcContent[[input$dataset]][["molPharmData"]][["exp"]]
		comparisonData <- comparisonData[, names(responseVec)]
		
		#pcGeneSets <- stringr::str_split(stringr::str_trim(input$pcGeneSets), pattern = "\\s+")[[1]]
		pcGeneSets <- input$pcGeneSets
		if (!("All Genes" %in% pcGeneSets)){
			pcGenes <- sort(unique(c(geneSetPathwayAnalysis::geneSets[pcGeneSets], recursive = TRUE)))
			pcGenes <- paste0("exp", pcGenes)
			pcGenes <- intersect(pcGenes, rownames(comparisonData))
			comparisonData <- comparisonData[pcGenes, ]
		}
		
		# ----[enable progress bar]--------------------------------------------------
		progress <- shiny::Progress$new()
		progress$set(message = "Computing Pattern Comparison Results: ", value = 0)
		# Close the progress when this reactive exits (even if there's an error).
		on.exit(progress$close())
		N <- nrow(comparisonData)
		updateProgress <- function(detail = NULL) {
			progress$inc(amount = 1/N, detail = detail)
		}
		# ---------------------------------------------------------------------------
		pcResults <- rcellminer::parCorPatternComparison(x = responseVec,
																										 Y = comparisonData,
																										 Z = currentPredictorData,
																										 updateProgress = updateProgress)
	
		return(pcResults)
	})
	
	#--------------------------------------------------------------------------------------
	
	if(require(rCharts)) {
		#----[Show 2D Actual vs. Predicted Response Scatter Plot in 'Plot' Tab]----------------
		output$plot <- renderChart({
			responseData <- rmResponseData()
			rmAlgoResults <- algoResults()
			predResponseData <- list()
			predResponseData$name <- paste0("predicted_", responseData$name)
			predResponseData$data <- rmAlgoResults$predictedResponse
			predResponseData$plotLabel <- predResponseData$name
			predResponseData$uniqName  <- predResponseData$name
			predResponseData$dataSource <- responseData$dataSource
			
			h1 <- makePlot(xData = predResponseData, yData = responseData, showColor = TRUE,
							showColorTissues = "all", dataSource = input$dataset,
							selectedTissuesOnly = FALSE, srcContent = srcContentReactive(), 
							dom="rm-plot", showPValue = FALSE)
		})
	
		#----[Show 2D Actual vs. CV-Predicted Response Scatter Plot in 'Cross-Validation' Tab]-
		output$cvPlot <- renderChart({
			responseData <- rmResponseData()
			rmAlgoResults <- algoResults()
			cvPredResponseData <- list()
			cvPredResponseData$name <- paste0("cv_predicted_", responseData$name)
			cvPredResponseData$data <- rmAlgoResults$cvPredictedResponse
			cvPredResponseData$plotLabel <- cvPredResponseData$name
			cvPredResponseData$uniqName  <- cvPredResponseData$name
			cvPredResponseData$dataSource <- responseData$dataSource
			
			h1 <- makePlot(xData = cvPredResponseData, yData = responseData, showColor = TRUE,
										 showColorTissues = "all", dataSource = input$dataset,
										 selectedTissuesOnly = FALSE, srcContent = srcContentReactive(), 
										 dom="rm-cvPlot", showPValue = FALSE)
		})
	}
	
	#----[Show Input + Predicted Response Data in 'Data' Tab]-------------------------------
	output$data <- DT::renderDataTable({
		dat <- inputData()
		rmAlgoResults <- algoResults()
		
		dat <- cbind(predicted_response = signif(rmAlgoResults$predictedResponse, 3), dat[, -1])
		if (length(rmAlgoResults$cvPredictedResponse) > 0){
			dat <- cbind(cv_predicted_response = signif(rmAlgoResults$cvPredictedResponse, 3), dat)
		}
		dat <- cbind(CellLine = rownames(dat), dat)
		
		DT::datatable(dat, rownames=FALSE, colnames=colnames(dat), filter='top', 
									style='bootstrap', options=list(pageLength = nrow(dat)))
	})
	
	#----[Show Predictors and Response in 'Heatmap' Tab]------------------------------------
	output$heatmap <- renderD3heatmap({
		dataTab <- inputData()
		dataMatrix <- as.matrix(t(dataTab[, -1]))
		numHiLoCols <- min(input$numHiLoResponseLines, floor(ncol(dataMatrix)/2))
		# Order columns by decreasing response values.
		dataMatrix <- dataMatrix[, order(dataMatrix[1, , drop = TRUE], decreasing = TRUE)]
		# Extract highest/lowest responder columns.
		colIndexSet <- c(1:numHiLoCols, (ncol(dataMatrix) - numHiLoCols + 1):ncol(dataMatrix))
		dataMatrix <- dataMatrix[, colIndexSet]
		# Remove data source identifier in predictor names.
		rownames(dataMatrix) <- vapply(rownames(dataMatrix), function(x) { 
			stringr::str_split(x, "_")[[1]][1] }, character(1))
		
		xAxisFontSize <- "6pt"
		# Notes:
		# Investigate 'cellnote' argument, which may allow use one matrix with suitably
		# scaled data for heatmap (color scaling), with original values displayed in 
		# tooltips.
		d3heatmap::d3heatmap(dataMatrix, dendrogram = "none", colors = "RdYlGn",
												 xaxis_font_size = xAxisFontSize, 
												 #width = (ncol(dataMatrix)*4),
												 xaxis_height = 200, yaxis_width = 200)
	})
	
	#----[Show Technical Details in 'Technical Details' Tab]--------------------------------
	output$techDetails <- renderPrint({
		rmAlgoResults <- algoResults()
		summary(rmAlgoResults$techDetails)
	})
	
	output$patternCompResults <- DT::renderDataTable({
		pcResults <- parCorPatternCompResults()
		
		DT::datatable(pcResults, rownames=FALSE, colnames=colnames(pcResults), filter='top', 
									style='bootstrap', 
									options=list(lengthMenu = c(10, 25, 50, 100), pageLength = 10))
	})
	
	#----[Organize Above Tabs for Display]--------------------------------------------------
	output$tabsetPanel = renderUI({
		maxNumHiLoResponseLines <- 100
		
		# TABS: Plot, Cross-Validation, Data, Heatmap, Technical Details
		ns <- session$ns
		dataTabPanel <- tabPanel("Data", DT::dataTableOutput(ns("data")))
		heatmapTabPanel <- tabPanel("Heatmap", 
																sliderInput(ns("numHiLoResponseLines"), 
																						"Number of High/Low Response Lines to Display:", 
																						min=1, max=maxNumHiLoResponseLines, 
																						value=20, width = "50%"),
																d3heatmapOutput(ns("heatmap")))
		techDetailsTabPanel <- tabPanel("Technical Details", verbatimTextOutput(ns("techDetails")))
		patternCompTabPanel <- tabPanel("Partial Correlation", 
																		selectInput(ns("pcGeneSets"), "Select Gene Sets",
																								choices  = c(names(geneSetPathwayAnalysis::geneSets), "All Genes"),
																								selected = "All Gene Sets",
																								multiple=TRUE),
																		actionButton(ns("computeParCors"), "Run"),
																		tags$hr(),
																		DT::dataTableOutput(ns("patternCompResults")))
		
		if (require(rCharts)){
			plotTabPanel   <- tabPanel("Plot", showOutput(ns("plot"), "highcharts"))
			cvPlotTabPanel <- tabPanel("Cross-Validation", showOutput(ns("cvPlot"), "highcharts"))
			tabsetPanel(type = "tabs", heatmapTabPanel, dataTabPanel, plotTabPanel, 
									cvPlotTabPanel, techDetailsTabPanel, patternCompTabPanel)
		} else{
			tabsetPanel(type = "tabs", heatmapTabPanel, dataTabPanel, techDetailsTabPanel, 
									patternCompTabPanel)	
		}
	})
	
	#-----[Support for Reactive UI Elements]-----------------------------------------------
	output$responseDataTypeUi <- renderUI({
		ns <- session$ns
		srcContent <- srcContentReactive()
		selectInput(ns("responseDataType"), "Response Data Type",
								choices  = srcContent[[input$dataset]][["featurePrefixes"]],
								selected = srcContent[[input$dataset]][["defaultFeatureY"]])
	})
	
	output$predDataTypesUi <- renderUI({
		ns <- session$ns
		srcContent <- srcContentReactive()
		selectInput(ns("predDataTypes"), "Predictor Data Types",
								choices  = srcContent[[input$dataset]][["featurePrefixes"]],
								selected = srcContent[[input$dataset]][["defaultFeatureX"]],
								multiple=TRUE)
	})
	#--------------------------------------------------------------------------------------
	
	#************************************************************************************************
}