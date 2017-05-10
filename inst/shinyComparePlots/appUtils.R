source("regressionModels.R")
#--------------------------------------------------------------------------------------------------
# Helper functions.
#--------------------------------------------------------------------------------------------------
getMatchedIds <- function(prefix, id, dataSource, srcContent){
	dat <- srcContent[[dataSource]][["molPharmData"]][[prefix]]
	idSet <- rcellminer::removeMolDataType(rownames(dat))
	
	if(id %in% idSet){
		return(id)
	}
	
	matchedIds <- NULL
	
	i <- which(toupper(id) == toupper(idSet))
	if (length(i) == 1){
		# Straightforward case-insensitive match.
		matchedIds <- unname(idSet[i])
	} else{
		# For drugs: try to match synonyms to source-specific identifiers.
		if (require(rcellminerUtils) && isDrugActivityDataType(prefix)){
			matchedIds <- rcellminerUtils::getDbDrugIds(drugName = id, dbName = dataSource)
			matchedIds <- intersect(matchedIds, idSet)
		}
	}
	
	return(matchedIds)
}

validateEntry <- function(prefix, id, dataSource, srcContent) {
	molPharmData <- srcContent[[dataSource]][["molPharmData"]]
	
	if(paste0(prefix, id) %in% rownames(molPharmData[[prefix]])) {
		return(TRUE)
	}
	
	return(FALSE)
}

getFeatureData <- function(prefix, id, dataSource, srcContent) {
	molPharmData <- srcContent[[dataSource]][["molPharmData"]]
	
	name <- paste0(prefix, id)
	data <- as.numeric(molPharmData[[prefix]][name, ])
	names(data) <- names(molPharmData[[prefix]][name, ])
	
	results <- list(name=name, data=data)
	
	# e.g., expTOP1 with dataSource=nci60 becomes TOP1 (exp, nci60)
	results$plotLabel <- paste0(id, " (", prefix, ", ", dataSource, ")")
	
	# e.g., expTOP1 with dataSource=nci60 becomes expTOP1_nci60; needed for 
	# getPlotData() results (data.frame) with data for same feature from different sources.
	results$uniqName <- paste0(results$name, "_", dataSource)
	
	results$dataSource <- dataSource
	
	return(results)
}

getTissueTypeSamples <- function(tissueTypes, dataSource, srcContent) {
	matchedSamples <- c(lapply(tissueTypes, function(tissue){
		srcContent[[dataSource]]$tissueToSamplesMap[[tissue]]
	}), recursive=TRUE)
	return(unique(matchedSamples))
}

# Returns all tissue types associated with one or more samples in sampleSet.
getSampleSetTissueTypes <- function(sampleSet, dataSource, srcContent) {
	tissueToSamples <- srcContent[[dataSource]]$tissueToSamplesMap
	
	isMatchedType <- vapply(names(tissueToSamples), function(tissueType) {
		length(intersect(sampleSet, tissueToSamples[[tissueType]])) > 0
	}, logical(1))
	
	matchedTypes <- character(0)
	if (any(isMatchedType)) {
		matchedTypes <- sort(unique(names(tissueToSamples[isMatchedType])))
	}
	
	return(matchedTypes)
}

getPlotData <- function(xData, yData, showColor, showColorTissues, dataSource=NULL,
												srcContent){
	if (is.null(dataSource)){
		dataSource <- xData$dataSource
	}
	
	#-----[make sure x and y data cell lines are matched]----------------------------------
	if (xData$dataSource == yData$dataSource){
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
	
	# RESTRICT (to rows with no NAs in either column x or column y).
	notNaData <- (!is.na(df[, xData$uniqName])) & (!is.na(df[, yData$uniqName]))
	df <- df[notNaData, ]
	
	if (nrow(df) > 0){
		cellLineSet <- rownames(df)
		# ADD COLOR COLUMN --------------------------------------------------------------------
		if (showColor){
			# Are there any cell lines to highlight in red?
			highlightedLineSet <- character(0)
			if (length(showColorTissues) > 0){
				highlightedLineSet <- unname(intersect(cellLineSet,
					getTissueTypeSamples(showColorTissues, dataSource, srcContent)))
			}
			
			if (length(highlightedLineSet) > 0){
				colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
				names(colorsToUse) <- rownames(df)
				colorsToUse[highlightedLineSet] <- "rgba(255,0,0,0.7)" # red
			} else{
				sampleTissueTypes <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
				colorsToUse <- srcContent[[dataSource]]$tissueColorMap[sampleTissueTypes]
			}
		} else{
			colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
		}
		df$color <- colorsToUse
		
		# ADD ONCOTREE TISSUE TYPE COLUMNS ----------------------------------------------------
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
		
		if ("EMT" %in% colnames(srcContent[[dataSource]]$sampleData)) {
			df$EMT <- srcContent[[dataSource]]$sampleData[rownames(df), "EMT"]
		}
	}
	
	return(df)
}


makePlot <- function(xData, yData, showColor, showColorTissues, dataSource,
										 srcContent, dom="rCharts", showPValue = TRUE) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, srcContent)
	
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
	
	# corResults <-cor.test(df[,xData$uniqName], df[,yData$uniqName], use="pairwise.complete.obs")
	# title <- paste0(paste(yData$plotLabel, '~', xData$plotLabel),
	# 								', r=', round(corResults$estimate, 2),
	# 								' p=', signif(corResults$p.value, 2))
	
	# rbind(y, y) forces more precise pvalue computation.
	corResults <- crossCors(df[,xData$uniqName], 
													rbind(df[,yData$uniqName], df[,yData$uniqName]))
	title <- paste0(paste(yData$plotLabel, '~', xData$plotLabel),
									', r=', signif(corResults$cor[1], digits=3))
	if (showPValue){
		title <- paste0(title, ' p=', signif(corResults$pval[1], digits=3))
	}
	
	h1$title(text=title)
	
	xAxisMin <- min(xData$data, na.rm = TRUE) - 0.25
	xAxisMax <- max(xData$data, na.rm = TRUE) + 0.25
	
	yAxisMin <- min(yData$data, na.rm = TRUE) - 0.25
	yAxisMax <- max(yData$data, na.rm = TRUE) + 0.25
	
	h1$xAxis(title=list(enabled=TRUE, text=xData$plotLabel, style=list(fontSize="24px", fontWeight="bold")),
					 min=xAxisMin, max=xAxisMax, labels=list(style=list(fontSize="20px")))
	h1$yAxis(title=list(enabled=TRUE, text=yData$plotLabel, style=list(fontSize="24px", fontWeight="bold")),
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
	
	h1$chart(zoomType="xy", style=list(fontFamily="Helvetica Neue"))
	
	# Enable exporting
	h1$exporting(enabled=TRUE)
	
	# Set name
	h1$set(dom=dom)
	
	# Print chart
	return(h1)
}


makePlotStatic <- function(xData, yData, showColor, showColorTissues, dataSource, 
													 srcContent, xLimVals = NULL, yLimVals = NULL) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, srcContent)
	df$tooltip <- paste0(
		"Cell: ", df$name, "\n",
		"Tissue: ", df$PlotTissueType, "\n",
		"OncoTree1: ", df$OncoTree1, "\n",
		"OncoTree2: ", df$OncoTree2, "\n",
		xData$uniqName, ": ", round(df$x, 2), "\n",
		yData$uniqName, ": ", round(df$y, 2), "\n"
	)
	
	# colorTab <- loadNciColorSet(returnDf=TRUE)
	# tissueColorTab <- unique(colorTab[, c("tissues", "colors")])
	tooltipCol <- "tooltip"
	
	# Plot parameters 
	classCol <- "color"
	colorPalette <- df[, "color"]
	names(colorPalette) <- df[, classCol]
	
	# Merge data
	df[, classCol] <- as.factor(df[, classCol])
	
	p1 <- rcellminer::plotCellMiner2D(df, xCol="x", yCol="y", xLabel = xData$plotLabel, yLabel = yData$plotLabel,
												colorPalette=colorPalette, classCol=classCol, tooltipCol=tooltipCol,
												xLimVal = xLimVals, yLimVal = yLimVals)
	
	return(p1)
}

getLmEquationString <- function(predictorWts, orderByDecrAbsVal = TRUE, numSigDigits = 3){
	if (length(predictorWts) == 0){
		return("")
	}
	if (is.null(names(predictorWts))){
		names(predictorWts) <- paste0("predictor_", 1:length(predictorWts))
	}
	if (orderByDecrAbsVal){
		predictorWts <- predictorWts[order(abs(predictorWts), decreasing = TRUE)]
	}
	if ((length(predictorWts) > 0) && ("(Intercept)" %in% names(predictorWts))){
		i <- which(names(predictorWts) == "(Intercept)")
		predictorWts <- c(predictorWts[i], predictorWts[-i])
	}
	
	predictorWts <- signif(predictorWts, numSigDigits)
	if (names(predictorWts)[1] == "(Intercept)"){
		eqStr <- paste0("Y = ", predictorWts[1])
	} else{
		eqStr <- paste0("Y = (", predictorWts[1], "*", names(predictorWts)[1], ")")
	}
	
	if (length(predictorWts) > 1){
		for (predName in names(predictorWts[-1])){
			eqStr <- paste0(eqStr, " + (", predictorWts[predName], "*", predName, ")")
		}
	}
		
	return(eqStr)
}

#--------------------------------------------------------------------------------------------------