#' Title
#'
#' @param sampleData
#' @param typeLevelSeparator
#'
#' @return value
#'
#' @examples
#' NULL
getTissueToSamplesMap <- function(sampleData, typeLevelSeparator = ":"){
	stopifnot(all(!duplicated(sampleData$Name)))
	rownames(sampleData) <- sampleData$Name
	tissueToSamples <- list()
	ocLevels <- paste0("OncoTree", 1:4)
	
	for (sample in rownames(sampleData)){
		sampleOcTypes <- as.character(sampleData[sample, ocLevels])
		typeName <- sampleOcTypes[1]
		if (is.na(typeName)){
			next
		}
		
		tissueToSamples[[typeName]] <- c(tissueToSamples[[typeName]], sample)
		for (i in (2:4)){
			if (is.na(sampleOcTypes[i])){
				break
			}
			typeName <- paste0(typeName, typeLevelSeparator, sampleOcTypes[i])
			tissueToSamples[[typeName]] <- c(tissueToSamples[[typeName]], sample)
		}
	}
	
	# ----[test]------------------------------------------------------------------
	# 	stopifnot(identical(sort(unique(c(tissueToSamples, recursive = TRUE))),
	# 						sort(rownames(sampleData))))
	# 	for (typeName in names(tissueToSamples)){
	# 		ocTypes <- str_split(typeName, pattern = typeLevelSeparator)[[1]]
	#
	# 		for (sample in tissueToSamples[[typeName]]){
	# 			sampleOcTypes <- as.character(sampleData[sample, ocLevels])
	# 			stopifnot(identical(sampleOcTypes[1:length(ocTypes)], ocTypes))
	# 		}
	# 	}
	# ----------------------------------------------------------------------------
	
	return(tissueToSamples)
}


#' Title
#'
#' @param dataPkgName
#'
#' @return value
#'
#' @examples
#' NULL
loadSourceContent <- function(dataPkgName){
	if (!require(dataPkgName, character.only = TRUE)){
		stop(paste0("Package '", dataPkgName, "' is not available."))
	}
	srcEnv <- new.env()
	data("molData", package=dataPkgName, verbose=TRUE, envir=srcEnv)
	data("drugData", package=dataPkgName, verbose=TRUE, envir=srcEnv)
	
	src <- list()
	src$molPharmData <- getAllFeatureData(srcEnv$molData)
	src$molPharmData[["act"]] <- exprs(getAct(srcEnv$drugData))
	
	for (featureType in names(src$molPharmData)){
		rownames(src$molPharmData[[featureType]]) <-
			paste0(featureType, rownames(src$molPharmData[[featureType]]))
	}
	
	# TO DO: Update to obtain this information from featureData(getAct(srcEnv$drugData))
	src$drugInfo <- data.frame(ID = rownames(exprs(getAct(srcEnv$drugData))), stringsAsFactors = FALSE)
	src$drugInfo$NAME <- src$drugInfo$ID
	src$drugInfo$MOA <- character(nrow(src$drugInfo))
	
	stopifnot(identical(unname(removeMolDataType(rownames(src$molPharmData$act))),
											src$drugInfo$ID))
	rownames(src$drugInfo) <- rownames(src$molPharmData$act)
	
	src$sampleData <- getSampleData(srcEnv$molData)
	rownames(src$sampleData) <- src$sampleData$Name
	
	# TO DO: Check whether spaces in tissue sample names creates any problems.
	src$tissueToSamplesMap <- getTissueToSamplesMap(src$sampleData)
	
	# TO DO: Properly define color map.
	src$tissueColorMap <- rep("rgba(0,0,255,0.5)", length(src$tissueToSamplesMap))
	names(src$tissueColorMap) <- names(src$tissueToSamplesMap)
	
	return(src)
}

#--------------------------------------------------------------------------------------------------
# Helper functions.
#--------------------------------------------------------------------------------------------------

#' Title
#'
#' @param prefix
#' @param id
#' @param dataSource
#' @param srcContent
#'
#' @return value
#'
#' @examples
#' NULL
validateEntry <- function(prefix, id, dataSource, srcContent) {
	molPharmData <- srcContent[[dataSource]][["molPharmData"]]
	
	if(paste0(prefix, id) %in% rownames(molPharmData[[prefix]])) {
		return(TRUE)
	}
	
	return(FALSE)
}

#' Title
#'
#' @param prefix
#' @param id
#' @param dataSource
#' @param srcContent
#'
#' @return value
#'
#' @examples
#' NULL
getFeatureData <- function(prefix, id, dataSource, srcContent) {
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

#' Title
#'
#' @param xData
#' @param yData
#' @param showColor
#' @param showColorTissues
#' @param dataSource
#' @param srcContent
#'
#' @return value
#'
#' @examples
#' NULL
getPlotData <- function(xData, yData, showColor, showColorTissues, dataSource=NULL, srcContent){
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
	if(showColor) {
		if("all" %in% showColorTissues) {
			sampleTissueTypes <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
			colorsToUse <- srcContent[[dataSource]]$tissueColorMap[sampleTissueTypes]
		} else {
			colorsToUse <- rep("rgba(0,0,255,0.5)", nrow(df)) #blue
			names(colorsToUse) <- rownames(df)
			
			for (tissueType in showColorTissues){
				matchedSamples <- intersect(srcContent[[dataSource]]$tissueToSamplesMap[[tissueType]],
																		rownames(df))
				colorsToUse[matchedSamples] <- "rgba(255,0,0,0.7)" # red
			}
		}
	} else{
		colorsToUse <- rep("rgba(0,0,255,0.5)", nrow(df)) #blue
	}
	
	df$color <- colorsToUse
	
	# Restrict to rows with no NAs in either column x or column y.
	notNaData <- (!is.na(df[, xData$uniqName])) & (!is.na(df[, yData$uniqName]))
	df <- df[notNaData, ]
	
	return(df)
}

#' Title
#'
#' @param xData
#' @param yData
#' @param showColor
#' @param showColorTissues
#' @param dataSource
#' @param srcContent
#'
#' @return value
#'
#' @examples
#' NULL
makePlot <- function(xData, yData, showColor, showColorTissues, dataSource, srcContent) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, srcContent)
	
	# Scatter plot
	h1 <- rCharts::Highcharts$new()
	
	# Divide the dataset, split by category and put into list() format
	# From: http://rcharts.io/viewer/?5735146#.VF6NS4W1Fy4
	series <- lapply(split(df, df$tissues), function(x) {
		res <- lapply(split(x, rownames(x)), as.list)
		names(res) <- NULL
		return(res)
	})
	
	invisible(sapply(series, function(x) {
		h1$series(data=x, type="scatter", name=x[[1]]$tissues)
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
	h1$xAxis(title=list(enabled=TRUE, text=xData$plotLabel, style=list(fontSize="20px", fontWeight="bold")),
					 labels=list(style=list(fontSize="20px")))
	h1$yAxis(title=list(enabled=TRUE, text=yData$plotLabel, style=list(fontSize="20px", fontWeight="bold")),
					 labels=list(style=list(fontSize="20px")))
	
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


#' Title
#'
#' @param xData
#' @param yData
#' @param showColor
#' @param showColorTissues
#' @param dataSource
#' @param srcContent
#'
#' @return value
#'
#' @examples
#' NULL
makePlotStatic <- function(xData, yData, showColor, showColorTissues, dataSource, srcContent) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, srcContent)
	
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