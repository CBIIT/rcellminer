library(shiny)
library(rcellminer)
library(rcellminerData)
library(jsonlite)

# Note: The jsonlite package has a validate() function whose name clashes with a 
# validate() function provided by the shiny package. 
# shiny::validate() must be used in the following code.

#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND CONSTRUCT APP (MOLECULAR + DRUG) DATA.
#--------------------------------------------------------------------------------------------------
source(system.file("shinyComparePlots", "functions.R", package="rcellminer"))
config <- jsonlite::fromJSON(system.file("shinyComparePlots", "config.json", package="rcellminer"))

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

#--------------------------------------------------------------------------------------------------
# Helper functions.
#--------------------------------------------------------------------------------------------------

validateEntry <- function(prefix, id, dataSource) {
  molPharmData <- srcContent[[dataSource]][["molPharmData"]]
  
  if(paste0(prefix, id) %in% rownames(molPharmData[[prefix]])) {
    return(TRUE)
  }

  return(FALSE)
}

getFeatureData <- function(prefix, id, dataSource) {
	molPharmData <- srcContent[[dataSource]][["molPharmData"]]
	
	entry <- paste0(prefix, id)
	data <- as.numeric(molPharmData[[prefix]][entry, ])
	names(data) <- names(molPharmData[[prefix]][entry, ])
		
	results <- list(entry=entry, data=data)
	
	# e.g., expTOP1 with dataSource=nci60 becomes TOP1 (exp, nci60)
	results$plotLabel <- paste0(id, " (", prefix, ", ", dataSource, ")")
	
	# e.g., expTOP1 with dataSource=nci60 becomes expTOP1_nci60; needed for 
	# getPlotData() results (data.frame) with data for same feature from different sources.
	results$uniqName <- paste0(results$entry, "_", dataSource)
	
	results$dataSource <- dataSource
	
	return(results)
}

getTissueTypeSamples <- function(tissueTypes, dataSource, srcContent){
	matchedSamples <- c(lapply(tissueTypes, function(tissue){
		srcContent[[dataSource]]$tissueToSamplesMap[[tissue]]
	}), recursive=TRUE)
	return(unique(matchedSamples))
}


getPlotData <- function(xData, yData, showColor, showColorTissues, dataSource=NULL,
												selectedTissuesOnly){
	if (is.null(dataSource)){
		dataSource <- xData$dataSource
	}
	
	#-----[make sure x and y data cell lines are matched]----------------------------------
	if (xData$dataSource != yData$dataSource){
		if (require(rcellminerUtils)){
			matchedLinesTab <- getMatchedCellLines(c(xData$dataSource, yData$dataSource))
			xData$data <- xData$data[matchedLinesTab[, 1]]
			yData$data <- yData$data[matchedLinesTab[, 2]]
		} else{
			stop("Install rcellminerUtils package to find matched cell lines between different data sources.")
		}
	} else{
		stopifnot(identical(names(xData$data), names(yData$data)))
	}
	#--------------------------------------------------------------------------------------
	
	df <- data.frame(x=names(xData$data), y=xData$data, z=yData$data, stringsAsFactors = FALSE)
	rownames(df) <- df$x
	colnames(df) <- c("Cell Line", xData$uniqName, yData$uniqName)
	
	# HighCharts series name
  df$tissues <- srcContent[[dataSource]]$sampleData[rownames(df), "TissueType"]
  
	# HighCharts point name
	df$name <- srcContent[[dataSource]]$sampleData[rownames(df), "Name"]
  
	# HighCharts point x, y for scatter plot
	df$x <- df[,xData$uniqName]
	df$y <- df[,yData$uniqName]
  
  # NOTE: making assumption that tissue type sets are disjoint, which may not hold
  # once hierarchy of tissue types is introduced (OK, i.e., disjoint at OncoTree1 level).
	if (showColor) {
	  if (("all" %in% showColorTissues) || selectedTissuesOnly) {
	    sampleTissueTypes <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
	    colorsToUse <- srcContent[[dataSource]]$tissueColorMap[sampleTissueTypes]
	  } else {
	    colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
	    names(colorsToUse) <- rownames(df)
	    
	    matchedSamples <- getTissueTypeSamples(showColorTissues, dataSource, srcContent)
	    matchedSamples <- intersect(matchedSamples, rownames(df))
	    colorsToUse[matchedSamples] <- "rgba(255,0,0,0.7)" # red
	  }
	} else{
	  colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
	}
	
	df$color <- colorsToUse
	
	# Restrict to rows with no NAs in either column x or column y.
	notNaData <- (!is.na(df[, xData$uniqName])) & (!is.na(df[, yData$uniqName]))
	df <- df[notNaData, ]
	
	if (selectedTissuesOnly){
		if((length(showColorTissues) > 0) && (!("all" %in% showColorTissues))){
			matchedSamples <- getTissueTypeSamples(showColorTissues, dataSource, srcContent)
			df <- df[intersect(matchedSamples, rownames(df)), ]
		}
	}
	
	# Add OncoTree tissue type information, and guaranteed non-NA tissue type for 
	# use in plotting code, etc.
	df$OncoTree1 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
	df$OncoTree2 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree2"]
	df$OncoTree3 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree3"]
	df$OncoTree4 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree4"]
	df$PlotTissueType <- ifelse(is.na(df$tissues),
															paste0(df$OncoTree1, ifelse(is.na(df$OncoTree2), "", 
																													paste0(":", df$OncoTree2))), 
															df$tissues)
	if (any(is.na(df$PlotTissueType))){
		df$PlotTissueType[which(is.na(df$PlotTissueType))] <- "TISSUE_TYPE_NA"
	}
  	
	return(df)
}


makePlot <- function(xData, yData, showColor, showColorTissues, dataSource, selectedTissuesOnly) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, selectedTissuesOnly)

	# Scatter plot
	h1 <- rCharts::Highcharts$new()

	# Divide the dataset, split by category and put into list() format
	# From: http://rcharts.io/viewer/?5735146#.VF6NS4W1Fy4
	series <- lapply(split(df, df$PlotTissueType), function(x) {
		res <- lapply(split(x, rownames(x)), as.list)
		names(res) <- NULL
		return(res)
	})

	invisible(sapply(series, function(x) {
		h1$series(data=x, type="scatter", name=x[[1]]$PlotTissueType)
	}
	))

	# Regression Line
	fit <- lm('y~x', data=df)
	x1 <- min(df[[xData$uniqName]])
	y1 <- fit$coefficients[[2]]*x1 + fit$coefficients[[1]]
	x2 <- max(df[[xData$uniqName]])
	y2 <- fit$coefficients[[2]]*x2 + fit$coefficients[[1]]

	h1$series(data=list(c(x1,y1), c(x2,y2)), type="line", color="#FF0000",
						marker=list(enabled=FALSE), enableMouseTracking=FALSE)

	corResults <-cor.test(df[,xData$uniqName], df[,yData$uniqName], use="pairwise.complete.obs")

	title <- paste0(paste(yData$plotLabel, '~', xData$plotLabel),
									', r=', round(corResults$estimate, 2),
									' p=', signif(corResults$p.value, 2))

	h1$title(text=title)
	
	xAxisMin <- min(xData$data, na.rm = TRUE) - 0.25
	xAxisMax <- max(xData$data, na.rm = TRUE) + 0.25
	
	yAxisMin <- min(yData$data, na.rm = TRUE) - 0.25
	yAxisMax <- max(yData$data, na.rm = TRUE) + 0.25
		
	h1$xAxis(title=list(enabled=TRUE, text=xData$plotLabel, style=list(fontSize="20px", fontWeight="bold")),
					 min=xAxisMin, max=xAxisMax, labels=list(style=list(fontSize="20px")))
	h1$yAxis(title=list(enabled=TRUE, text=yData$plotLabel, style=list(fontSize="20px", fontWeight="bold")),
					 min=yAxisMin, max=yAxisMax, labels=list(style=list(fontSize="20px")))

	h1$legend(enabled=FALSE)

	# Force circle markers, set default size, hover color (otherwise color unpredictable)
	h1$plotOptions(series=list(animation=50),
		scatter=list(marker=list(symbol='circle', radius=6,
				 	 states=list(hover=list(fillColor='white')))))

	tooltipFormat <- paste0("#! function() { return 'Cell: ' + this.point.name +
		'<br/>Tissue: ' + this.series.name +
		'<br/>", xData$uniqName, ": ' + Math.round(this.x * 100) / 100 +
		'<br/>", yData$uniqName, ": ' + Math.round(this.y * 100) / 100; } !#")

	h1$tooltip(backgroundColor="rgba(255,255,255,1)", formatter=tooltipFormat)

	h1$chart(zoomType="xy")

	# Enable exporting
	h1$exporting(enabled=TRUE)

	# Set name
	h1$set(dom="rCharts")

	# Print chart
	return(h1)
}


makePlotStatic <- function(xData, yData, showColor, showColorTissues, dataSource, selectedTissuesOnly=FALSE) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, selectedTissuesOnly=FALSE)
	
	corResults <-cor.test(df[,xData$uniqName], df[,yData$uniqName], use="pairwise.complete.obs")
	
	title <- paste0(paste(yData$plotLabel, '~', xData$plotLabel),
									', r=', round(corResults$estimate, 2),
									' p=', signif(corResults$p.value, 2))
	
	plot(df[, xData$uniqName], df[, yData$uniqName], xlab=xData$plotLabel, ylab=yData$plotLabel, 
			 col=df[,"color"], pch=16, main=title)
	formula <- as.formula(paste(yData$uniqName, "~", xData$uniqName))
	abline(lm(formula, df), col="red")
}

#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
shinyServer(function(input, output, session) {
  
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
				shiny::validate(
		            need(validateEntry(input$xPrefix, input$xId, input$xDataset), 
                     paste("ERROR:", paste0(input$xPrefix, input$xId), "not found.")),
		            need(validateEntry(input$yPrefix, input$yId, input$yDataset), 
                     paste("ERROR:", paste0(input$yPrefix, input$yId), "not found."))
		    )
		
				xData <- getFeatureData(input$xPrefix, input$xId, input$xDataset)
				yData <- getFeatureData(input$yPrefix, input$yId, input$yDataset)
				
				h1 <- makePlot(xData, yData, input$showColor, input$showColorTissues, input$xDataset, 
											 input$selectedTissuesOnly)
				#return(h1)
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
			shiny::validate(
				need(validateEntry(input$xPrefix, input$xId, input$xDataset), 
						 paste("ERROR:", paste0(input$xPrefix, input$xId), "not found.")),
				need(validateEntry(input$yPrefix, input$yId, input$yDataset), 
						 paste("ERROR:", paste0(input$yPrefix, input$yId), "not found."))
			)
			
			xData <- getFeatureData(input$xPrefix, input$xId, input$xDataset)
			yData <- getFeatureData(input$yPrefix, input$yId, input$yDataset)
			
			makePlotStatic(xData, yData, input$showColor, input$showColorTissues, input$xDataset,
										 input$selectedTissuesOnly)
		})
		#--------------------------------------------------------------------------------------
    
    #----[Render Data Table in 'Download Data' Tab]----------------------------------------
    # Generate an HTML table view of the data
    output$table <- renderDataTable({
    	if (input$selectedTissuesOnly){
    		shiny::validate(need(length(input$showColorTissues) > 0, "Please select tissue types."))
    	}
    	if (!require(rcellminerUtils)){
    		shiny::validate(need(input$xDataset == input$yDataset, "ERROR: x and y axis data sets must be the same."))
    	}
      shiny::validate(need(validateEntry(input$xPrefix, input$xId, input$xDataset), "ERROR: x-Axis Entry Not Found."))
      shiny::validate(need(validateEntry(input$yPrefix, input$yId, input$yDataset), "ERROR: y-Axis Entry Not Found."))

    	xData <- getFeatureData(input$xPrefix, input$xId, input$xDataset)
    	yData <- getFeatureData(input$yPrefix, input$yId, input$yDataset)
			
    	# Column selection below is to restrict to cell line, x, y features, 
    	# and tissue type information (source-provided + OncoTree).
    	dlDataTab <- getPlotData(xData, yData, input$showColor, input$showColorTissues, 
    													 input$xDataset, input$selectedTissuesOnly)
    	dlDataTabCols <- c(colnames(dlDataTab)[1:4], paste0("OncoTree", 1:4))
    	dlDataTab <- dlDataTab[, dlDataTabCols]
    }, options = list(paging=FALSE))
		#--------------------------------------------------------------------------------------
    
    #----[Render Data Table in 'Search IDs' Tab]-------------------------------------------
    # Generate an HTML table view of the data
    output$ids <- renderDataTable({
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

        results

    }, options = list(paging=TRUE, pageLength=10))
		#--------------------------------------------------------------------------------------
    
		#----[Render Data Table in 'Compare Patterns' Tab]-------------------------------------
		output$patternComparison <- renderDataTable({
			if (input$selectedTissuesOnly){
				shiny::validate(need(length(input$showColorTissues) > 0, "Please select tissue types."))
			}
			if (!require(rcellminerUtils)){
				shiny::validate(need(input$xDataset == input$yDataset, "ERROR: x and y axis data sets must be the same."))
			}
		  shiny::validate(need(validateEntry(input$xPrefix, input$xId, input$xDataset), "ERROR: x-Axis Entry Not Found."))
		  
		  dat <- getFeatureData(input$xPrefix, input$xId, input$xDataset)
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
		    
		    # Reorder columns
		    results <- results[, c("ids", "NAME", "COR", "PVAL")]
		    colnames(results) <- c("ID", "Name", "Correlation", "P-Value")
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
		  
		  results
		}, options=list(pageLength=10))
		#--------------------------------------------------------------------------------------
    
		output$log <- renderText({
				paste(names(input), collapse=" ")
				query <- parseQueryString(session$clientData$url_search)
				sapply(names(input), function(x) { query[[x]] <- input[[x]] })
			})

		output$genUrl <- renderText({
			#query <- parseQueryString(session$clientData$url_search)
			query <- list()

			# Get current input values and combine with
			tmp <- sapply(names(input), function(x) { query[[x]] <- input[[x]] })
			query <- c(query, tmp)

			paramStr <- paste(sapply(names(query), function(x) { paste0(x, "=", query[[x]]) }), collapse="&")

			urlStr <- paste(session$clientData$url_protocol, "//",
											session$clientData$url_hostname, ":",
											session$clientData$url_port,
											session$clientData$url_pathname,
											"?", paramStr,
											sep="")

			paste0(a(paste("Shareable Link (Right-click then 'Copy' to share)"), href=urlStr), hr())

			#paste("Shareable Link:", urlStr)
		})

    #**********************************************************************************************
		output$tabsetPanel = renderUI({		
			#verbatimTextOutput("log") can be used for debugging
			#tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),
			
			tab1 <- tabPanel("Download Data", 
                       downloadLink("downloadData", "Download Data as Tab-Delimited File"), 
                       dataTableOutput("table"))
			tab2 <- tabPanel("Search IDs", 
                       includeMarkdown("www/files/help.md"), 
                       dataTableOutput("ids"))
			tab3 <- tabPanel("Compare Patterns", 
											 includeMarkdown("www/files/help.md"), 
                       selectInput("patternComparisonType", "Pattern Comparison to x-Axis Entry", 
                                   choices=c("Molecular Data"="molData", "Drug Data"="drug"), selected="molData"), 
                       dataTableOutput("patternComparison")) 
			tab4 <- tabPanel("About", 
			                 includeMarkdown("www/files/about.md"))
			
			if(input$hasRCharts == "TRUE") {
				tabsetPanel(type="tabs", 
										#tabPanel("Plot Data", htmlOutput("genUrl"), showOutput("rCharts", "highcharts")),
										tabPanel("Plot Data", showOutput("rCharts", "highcharts")),
										tab1, tab2, tab3, tab4					
				)				
			} else {
				tabsetPanel(type="tabs", 
										#tabPanel("Plot Data", htmlOutput("genUrl"), plotOutput("rChartsAlternative", width=600, height=600)),
										tabPanel("Plot Data", plotOutput("rChartsAlternative", width=600, height=600)),
										tab1, tab2, tab3, tab4
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
      	xData <- getFeatureData(input$xPrefix, input$xId, input$xDataset)
      	yData <- getFeatureData(input$yPrefix, input$yId, input$yDataset)
				
      	df <- getPlotData(xData, yData, input$showColor, input$showColorTissues, 
      										input$xDataset, input$selectedTissuesOnly)
      	
      	# Column selection below is to restrict to cell line, x, y features, 
      	# and tissue type information (source-provided + OncoTree).
      	dfCols <- c(colnames(df)[1:4], paste0("OncoTree", 1:4))
      	df <- df[, dfCols]

      	write.table(df, file, quote=FALSE, row.names=FALSE, sep="\t")
      }
    )
  
  output$xPrefixUi <- renderUI({          
  	selectInput("xPrefix", "x-Axis Type", 
  							choices = srcContent[[input$xDataset]][["featurePrefixes"]], 
  							selected = srcContent[[input$xDataset]][["defaultFeatureX"]])
  })

  output$yPrefixUi <- renderUI({
  	selectInput("yPrefix", "y-Axis Type", 
  							choices = srcContent[[input$yDataset]][["featurePrefixes"]], 
  							selected = srcContent[[input$yDataset]][["defaultFeatureY"]])
  })

  output$showColorTissuesUi <- renderUI({
  	tissueTypes <- names(srcContent[[input$xDataset]][["tissueToSamplesMap"]])
  	selectInput("showColorTissues", "Color Specific Tissues?", 
  							choices=c("all", unique(tissueTypes)), multiple=TRUE, selected="all")
  })

})
#-----[end of shinyServer()]-----------------------------------------------------------------------
