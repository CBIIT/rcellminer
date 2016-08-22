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
	rmResponseData <- reactive({
		shiny::validate(need(validateEntry(input$responseDataType, 
																			 input$responseId, input$dataset,
																			 srcContent = srcContentReactive()),
												 paste("ERROR:", paste0(input$responseDataType, input$responseId), 
												 			"not found.")))
		yData <- getFeatureData(input$responseDataType, input$responseId, input$dataset, 
														srcContent = srcContentReactive())
		return(yData)
	})
	
	# Returns a data frame with data suitable for predictive modeling (updated according
	# to current user selections). The first column indicates the cell line, with
	# subsequent columns containing response and predictor variables.
	inputData <- reactive({
		yData <- rmResponseData()
		dataTab <- data.frame(CellLine = names(yData$data), stringsAsFactors = FALSE)
		rownames(dataTab) <- dataTab$CellLine
		dataTab[, yData$uniqName] <- yData$data
		
		predIds <- stringr::str_split(stringr::str_trim(input$predIds),pattern = "\\s+")[[1]] 
		for (id in predIds){
			for (dataType in input$predDataTypes){
				if (validateEntry(dataType, id, input$dataset, srcContentReactive())){
					xData <- getFeatureData(dataType, id, input$dataset, srcContentReactive())
					dataTab[, xData$uniqName] <- xData$data
				} else{
					warning(paste0(dataType, id), " not found.")
				}
			}
		}
		
		shiny::validate(need(ncol(dataTab) > 2,
												 paste("ERROR: No data for specified predictors.")))
		
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
			lmCvFit <- rcellminerElasticNet::getLmCvFit(X = as.matrix(lmData[, -1]), 
																									y = lmData[, 1, drop = TRUE], nFolds = 10, nRepeats = 1)
			# -----[assemble results]---------------------------------------------------
			rmAlgoResults$algorithm <- "Linear Regression"
			rmAlgoResults$predictorWts <- coef(lmFit)
			rmAlgoResults$predictedResponse <- predict(lmFit)
			rmAlgoResults$cvPredictedResponse <- lmCvFit$cvPred
			rmAlgoResults$techDetails <- lmFit
			# --------------------------------------------------------------------------
		} else{
			shiny::validate(FALSE, paste("ERROR: Algorithm not available."))
		}
		
		return(rmAlgoResults)
	})
	
	#--------------------------------------------------------------------------------------
	
	#----[Render 2D Plot in 'Plot' Tab]----------------------------------------------------
	if(require(rCharts)) {
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
							selectedTissuesOnly = FALSE, srcContent = srcContentReactive(), dom="rm-plot")
		})
	}
	
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
	
	output$heatmap <- renderD3heatmap({
		dataTab <- inputData()
		dataMatrix <- as.matrix(t(dataTab[, -1]))
		# Order columns by increasing response values.
		dataMatrix <- dataMatrix[, order(dataMatrix[1, , drop = TRUE])]
		# Remove data source identifier in predictor names.
		rownames(dataMatrix) <- vapply(rownames(dataMatrix), function(x) { 
			stringr::str_split(x, "_")[[1]][1] }, character(1))
		
		d3heatmap::d3heatmap(dataMatrix, dendrogram = "none", colors = "RdYlGn")
	})
	
	output$techDetails <- renderPrint({
		rmAlgoResults <- algoResults()
		summary(rmAlgoResults$techDetails)
	})
	
	output$tabsetPanel = renderUI({
		# TABS: Plot, Cross-Validation, Data, Heatmap, Technical Details
		ns <- session$ns
		dataTabPanel <- tabPanel("Data", DT::dataTableOutput(ns("data")))
		heatmapTabPanel <- tabPanel("Heatmap", d3heatmapOutput(ns("heatmap")))
		techDetailsTabPanel <- tabPanel("Technical Details", verbatimTextOutput(ns("techDetails")))
		
		if (require(rCharts)){
			plotTabPanel <- tabPanel("Plot", showOutput(ns("plot"), "highcharts"))
			tabsetPanel(type = "tabs", plotTabPanel,
									dataTabPanel, heatmapTabPanel, techDetailsTabPanel)
		} else{
			tabsetPanel(type = "tabs", dataTabPanel, heatmapTabPanel, techDetailsTabPanel)	
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