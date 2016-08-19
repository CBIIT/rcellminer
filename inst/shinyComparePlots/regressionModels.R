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
					 				selectInput(ns("rmDataset"), "Dataset", choices=dataSourceChoices, selected = "nci60"),
					 				uiOutput(ns("rmResponseDataTypeUi")),
					 				textInput(ns("rmResponseId"), "Response ID: (Case-Sensitive, e.g., 609699)", "609699"),
					 				uiOutput(ns("rmPredDataTypesUi")),
					 				textInput(ns("rmPredIds"), "Predictor IDS: (Case-Sensitive, e.g. SLFN11 JAG1)", "SLFN11 JAG1"),
					 				selectInput(ns("rmAlgorithm"), "Algorithm", 
					 										choices=c("Linear Regression", "Supervised Principal Components"), 
					 										selected = "Linear Regression")
					 			)
					 		),
					 		mainPanel(
					 			uiOutput(ns('rmTabsetPanel'))
					 		)
					 	)
					 )						 
	)
}
#-----[NavBar Tab: Regression Models (Server code)]------------------------------------------------
regressionModels <- function(input, output, session, srcContentReactive) {

	#----[Reactive Variables]--------------------------------------------------------------
	# Returns a data frame with data suitable for predictive modeling (updated according
	# to current user selections). The first column indicates the cell line, with
	# subsequent columns containing response and predictor variables.
	rmInputData <- reactive({
		shiny::validate(need(validateEntry(input$rmResponseDataType, 
																			 input$rmResponseId, input$rmDataset,
																			 srcContent = srcContentReactive()),
												 paste("ERROR:", paste0(input$rmResponseDataType, input$rmResponseId), 
												 			"not found.")))
		yData <- getFeatureData(input$rmResponseDataType, input$rmResponseId, input$rmDataset, 
														srcContent = srcContentReactive())
		dataTab <- data.frame(CellLine = names(yData$data), stringsAsFactors = FALSE)
		rownames(dataTab) <- dataTab$CellLine
		dataTab[, yData$uniqName] <- yData$data
		
		predIds <- stringr::str_split(stringr::str_trim(input$rmPredIds),pattern = "\\s+")[[1]] 
		for (id in predIds){
			for (dataType in input$rmPredDataTypes){
				if (validateEntry(dataType, id, input$rmDataset, srcContentReactive())){
					xData <- getFeatureData(dataType, id, input$rmDataset, srcContentReactive())
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
	rmAlgoResults <- reactive({
		dataTab <- rmInputData()
		algoResults <- list()
		
		# Note: refactor this function to gather and validate parameters, which are then passed
		# to a specialized implementation function that returns a standard format algoResults
		# list object.
		if (input$rmAlgorithm == "Linear Regression"){
			lmData <- dataTab[, -1] # First column has cell line names
			# Note: will have to handle issues caused by predictor names with spaces 
			# or other characters that cannot be used within a formula.
			lmFormula <- as.formula(paste0(colnames(lmData)[1], " ~ ."))
			lmFit <- lm(lmFormula, lmData)
			lmCvFit <- rcellminerElasticNet::getLmCvFit(X = as.matrix(lmData[, -1]), 
																									y = lmData[, 1, drop = TRUE], nFolds = 10, nRepeats = 1)
			# -----[assemble results]---------------------------------------------------
			algoResults$algorithm <- "Linear Regression"
			algoResults$predictorWts <- coef(lmFit)
			algoResults$predictedResponse <- predict(lmFit)
			algoResults$cvPredictedResponse <- lmCvFit$cvPred
			algoResults$techDetails <- lmFit
			# --------------------------------------------------------------------------
		} else{
			shiny::validate(FALSE, paste("ERROR: Algorithm not available."))
		}
		
		return(algoResults)
	})
	
	#--------------------------------------------------------------------------------------
	
	output$rmData <- DT::renderDataTable({
		dat <- rmInputData()
		cat("In DT::renderDataTable()")
		algoResults <- rmAlgoResults()
		
		dat <- cbind(predicted_response = signif(algoResults$predictedResponse, 3), dat[, -1])
		if (length(algoResults$cvPredictedResponse) > 0){
			dat <- cbind(cv_predicted_response = signif(algoResults$cvPredictedResponse, 3), dat)
		}
		dat <- cbind(CellLine = rownames(dat), dat)
		
		DT::datatable(dat, rownames=FALSE, colnames=colnames(dat), filter='top', 
									style='bootstrap', options=list(pageLength = nrow(dat)))
	})
	
	output$rmHeatmap <- renderD3heatmap({
		dataTab <- rmInputData()
		dataMatrix <- as.matrix(t(dataTab[, -1]))
		# Order columns by increasing response values.
		dataMatrix <- dataMatrix[, order(dataMatrix[1, , drop = TRUE])]
		# Remove data source identifier in predictor names.
		rownames(dataMatrix) <- vapply(rownames(dataMatrix), function(x) { 
			stringr::str_split(x, "_")[[1]][1] }, character(1))
		
		d3heatmap::d3heatmap(dataMatrix, dendrogram = "none", colors = "RdYlGn")
	})
	
	output$rmTechDetails <- renderPrint({
		algoResults <- rmAlgoResults()
		summary(algoResults$techDetails)
	})
	
	output$rmTabsetPanel = renderUI({
		# TABS: Plot, Cross-Validation, Data, Heatmap, Technical Details
		ns <- session$ns
		tabsetPanel(type = "tabs",
								tabPanel("Data", DT::dataTableOutput(ns("rmData"))),
								tabPanel("Heatmap", d3heatmapOutput(ns("rmHeatmap"))),
								tabPanel("Technical Details", verbatimTextOutput(ns("rmTechDetails"))))
	})
	
	#-----[Support for Reactive UI Elements]-----------------------------------------------
	output$rmResponseDataTypeUi <- renderUI({
		ns <- session$ns
		srcContent <- srcContentReactive()
		selectInput(ns("rmResponseDataType"), "Response Data Type",
								choices  = srcContent[[input$rmDataset]][["featurePrefixes"]],
								selected = srcContent[[input$rmDataset]][["defaultFeatureY"]])
	})
	
	output$rmPredDataTypesUi <- renderUI({
		ns <- session$ns
		srcContent <- srcContentReactive()
		selectInput(ns("rmPredDataTypes"), "Predictor Data Types",
								choices  = srcContent[[input$rmDataset]][["featurePrefixes"]],
								selected = srcContent[[input$rmDataset]][["defaultFeatureX"]],
								multiple=TRUE)
	})
	#--------------------------------------------------------------------------------------
	
	#************************************************************************************************
}