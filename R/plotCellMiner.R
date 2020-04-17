#' Description: Produces CellMiner-like plots in R 
#' 
#' @param drugAct a matrix of drug activity values (cell lines as columns, drug entries as rows)
#' @param molData a list of matricies a molecular 
#' @param plots a vector of characters denoting the plots to include and the order (e.g. c("mut", "drug", "cop"). 
#'   Currently, supported entries mutations (mut), drug activities (drug), copy number variations (cop)
#' @param nsc a string NSC ID that will be plotted when a "drug" entry appears in the plots vector  
#' @param gene a string HUGO gene symbol for which the "mut", "cop", or "exp" plots will be produced if in plots vector 
#' @param features a vector of strings that provide the full IDs for elements to be plotted (e.g. mutCDK4 for CDK4 mutations). 
#'   This overwrites the nsc and gene parameters, but is needed in advanced plots that involve data that involves one-to-many 
#'   relationships (e.g. many entries for a given gene in the exome data) and a gene symbol is ambiguous.
#' @param sub a vector of strings with sub-titles for each plot
#' @param xLimits a 2 number vector with the the minimum and maximum X-axis values (default: -3,3 for Z-scores, 0,1 for binary entries) 
#' @param xLabel a string for the default X-axis label
#' @param extraPlot a list containing title, label, and values (numeric vector of length 60); only one extra plot can be included
#' @param verbose a boolean to show debugging information
#' @return None
#' 
#' @author Augustin Luna <augustin AT mail.nih.gov>
#' 
#' @importFrom gplots barplot2
#' 
#' @examples
#' drugAct <- exprs(getAct(rcellminerData::drugData))
#' molData <- getMolDataMatrices()
#' plots <- c("mut", "drug", "cop", "xai", "pro")
#' plotCellMiner(drugAct, molData, plots=plots, nsc="94600", gene="CDK4", verbose=FALSE)
#' 
#' plots <- c("mut", "xai", "cop", "cop", "cop", "cop")
#' plotCellMiner(drugAct, molData, plots=plots, nsc="94600", gene=c("CDK4", "TP53", 
#'   "BRAF", "GAPDH"), verbose=FALSE)
#' 
#' plotCellMiner(drugAct, molData, plots=NULL, nsc=NULL, features=c("mutCDK4", 
#'   "xaiCDK4", "exochr1:101704532_G_T", "mdaIS_P53_MUT", "mirhsa-miR-22", "proTP53_26_GBL00064"), 
#'   verbose=FALSE)
#' 
#' @concept rcellminer
#' @export
#' 
# @importFrom graphics par layout lcm graphics
#' @importFrom utils str  
plotCellMiner <- function(drugAct, molData, plots, nsc=NULL, gene=NULL, 
													features=NULL, sub=NULL, xLimits=NULL, xLabel=NULL, 
													extraPlot=NULL, verbose=FALSE) {
	
	acceptablePlotTypes <- c("drug", "cop", "exp", "xai", "exo", "mut", "mir", "pro", "mda")
	
	if(!all(plots %in% acceptablePlotTypes)) {
		stop(paste('"plots" may only include:', paste(acceptablePlotTypes, collapse=" ")))
	}
	
	plotInfo <- loadCellminerPlotInfo(returnDf=TRUE)
	
	mainTitle <- NULL
	values <- NULL
	
	# Make sure the nsc is a character string otherwise NSCs may be treated as numeric indicies
	nsc <- as.character(nsc)
	
	if(is.null(features)) {
		copPrefix <- "cop"
		expPrefix <- "exp"
		xaiPrefix <- "xai"
		exoPrefix <- "exo"
		mutPrefix <- "mut"
		mirPrefix <- "mir"
		proPrefix <- "pro"
		
		copData <- molData[[copPrefix]]
		expData <- molData[[expPrefix]]
		xaiData <- molData[[xaiPrefix]]
		exoData <- molData[[exoPrefix]]
		mutData <- molData[[mutPrefix]]
		mirData <- molData[[mirPrefix]]
		proData <- molData[[proPrefix]]
		
		drugIdx <- 1
		copIdx <- 1
		expIdx <- 1
		xaiIdx <- 1
		exoIdx <- 1
		mutIdx <- 1
		mirIdx <- 1
		proIdx <- 1
		
		for(p in plots) {
			if(p == "drug") {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, nsc[drugIdx]))
				values <- cbind(values, as.numeric(drugAct[nsc[drugIdx],]))
				drugIdx <- drugIdx + 1
			}
			
			if(p == copPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[copIdx]))
				values <- cbind(values, as.numeric(copData[paste0(copPrefix, gene[copIdx]),]))
				copIdx <- copIdx + 1
			}	
			
			if(p == expPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[expIdx]))
				values <- cbind(values, as.numeric(expData[paste0(expPrefix, gene[expIdx]),]))
				expIdx <- expIdx + 1
			}
			
			if(p == xaiPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[xaiIdx]))
				values <- cbind(values, as.numeric(xaiData[paste0(xaiPrefix, gene[xaiIdx]),]))
				xaiIdx <- xaiIdx + 1
			}	
			
			if(p == exoPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[exoIdx]))
				values <- cbind(values, as.numeric(exoData[paste0(exoPrefix, gene[exoIdx]),]))
				exoIdx <- exoIdx + 1
			}
			
			if(p == mutPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[mutIdx]))
				values <- cbind(values, as.numeric(mutData[paste0(mutPrefix, gene[mutIdx]),]))
				mutIdx <- mutIdx + 1
			}
			
			if(p == mirPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[mirIdx]))
				values <- cbind(values, as.numeric(mirData[paste0(mirPrefix, gene[mirIdx]),]))
				mirIdx <- mirIdx + 1
			}	
			
			if(p == proPrefix) {
				plotPrefix <- plotInfo[plotInfo$featureType == p, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, gene[proIdx]))
				values <- cbind(values, as.numeric(proData[paste0(proPrefix, gene[proIdx]),]))
				proIdx <- proIdx + 1
			}	
		}
	} else {
		plots <- NULL
		
		for(fea in features) {
			prefix <- substr(fea, 1, 3)
			feature <- substr(fea, 4, nchar(fea))
			
			if(grepl("^\\d+$", fea)) {
				plotPrefix <- plotInfo[plotInfo$featureType == "drug", "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, fea))
				values <- cbind(values, as.numeric(drugAct[fea,]))
			} else {
				plotPrefix <- plotInfo[plotInfo$featureType == prefix, "titlePrefix"]
				mainTitle <- c(mainTitle, paste0(plotPrefix, feature))
				values <- cbind(values, as.numeric(molData[[prefix]][paste0(prefix, feature),]))	
			}
			
			plots <- c(plots, prefix)
		}
	}

	### Initiate Graphic  
	# NOTE: To make longer gene names fit use cex.axis 
	# in Heatmap section to scale down font or MAYBE increase the 
	# left margin c(BOTTOM, LEFT, TOP, RIGHT) in Heatmap and GI50
	
  # Save default for resetting
  opar <- par(no.readonly = TRUE)
  
  # Generate layout for plots and set x-axis limits
	figures <- length(plots) 
	figureCols <- 3
	m <- NULL
	limits <- NULL

	# Add data for extra plot
	if(!is.null(extraPlot)) {
		values <- cbind(values, extraPlot[["values"]])
		mainTitle <- c(mainTitle, extraPlot[["title"]])
		figures <- figures + 1
	}
		
	# Set subtitle
	if(!is.null(sub)) {
		subTitle <- sub
	} else {
		subTitle <- rep(NULL, figures)
	}

	for(i in 1:figures) {
		tmp <- NULL
		
		for(j in 1:figureCols) {
			tmp <- cbind(tmp, rep(i, 10))
		}

		m <- cbind(m, tmp) 
		
		# Set plot x-axis
		valuesTmp <- rev(values[,i])

		if(i <= length(plots)) {
			# Limits 
			curXMin <- plotInfo[plotInfo$featureType == plots[i], "xMin"]
			curXMax <- plotInfo[plotInfo$featureType == plots[i], "xMax"]
			
			if(verbose) {
				cat(paste(plots[i], collapse=" "), " Y ", curXMin, "\n")
				cat(paste(plots[i], collapse=" "), " X ", curXMax, "\n")
			}
	
			# Round to the nearest increment
			if(is.na(curXMin)) {
				curXMin <- floor(min(valuesTmp, na.rm=TRUE)*2)/2
			}
			
			if(is.na(curXMax)) {
				curXMax <- ceiling(max(valuesTmp, na.rm=TRUE)*2)/2
			}
					
			if(plots[i] %in% limits$prefix) {
				prevXMin <- limits[limits$prefix == plots[i], "xMin"]
				prevXMax <- limits[limits$prefix == plots[i], "xMax"]
				
				if(prevXMin > curXMin) {
					limits[limits$prefix == plots[i], "xMin"] <- curXMin
				}
				
				if(prevXMax < curXMax) {
					limits[limits$prefix == plots[i], "xMax"] <- curXMax
				}
			} else {
				limits <- rbind(limits, data.frame(prefix=plots[i], xMin=curXMin, xMax=curXMax, stringsAsFactors=FALSE))
			}
		}
	}
	layout(m, heights = lcm(12)) #widths = lcm(4), 
	#layout.show() 

	# Labels
	yLabel <- ""
	
	par(las=1, cex=0.5)
	# Margins: Bottom, Left, Top, Right
	par(oma=c(5,4.5,4,4))
	# Label Margins: 
	par(mgp=c(4, 1, 0))
	
	for(i in 1:figures) {
		## SET X LABEL
		if(!is.null(extraPlot) & i == (length(plots)+1)) {
			xLabel <- extraPlot[["label"]]	
		} else {
			xLabel <- plotInfo[plotInfo$featureType == plots[i], "label"]
			
			# Fix for MDA features
			if(plots[i] == "mda" & !is.null(features)) {
				xLabel <- plotInfo[plotInfo$featureType == features[i], "label"]
			} 
		}

		# Match CellMiner plots with BR lines at the top
		colorSet <- loadNciColorSet(returnDf=TRUE)
		colors <- rev(colorSet$colors)
		labels <- rev(colorSet$abbrCellLines) 
	  valuesTmp <- rev(values[,i])
	  
	  ## SET X-AXIS
	  if(!is.null(extraPlot) & i == figures) {
	  	xMin <- min(valuesTmp, na.rm=TRUE)
	  	xMax <- max(valuesTmp, na.rm=TRUE)
	  } else {
	  	xMin <- limits[limits$prefix == plots[i], "xMin"]
	  	xMax <- limits[limits$prefix == plots[i], "xMax"]	  
	  }
	  
	  #plotNci60Profile(values[,i], xMin, xMax, mainTitle[i], xLabel, yLabel, as.character(colors), labels)  
                     
	  if(verbose) {
	  	str(figures) 
	  	str(plots[i])
	  	str(labels) 
	  	str(values[,i]) 
	  	str(xLabel)
	  	str(xMin)
	  	str(xMax)
	  	str(limits)
	  }
	  
	  ## PLOT
		# Return the midpoints 
		mp <- barplot2(valuesTmp, axes=FALSE, axisnames=FALSE, xlim=c(xMin, xMax), 
                   main=mainTitle[i], sub=subTitle[i], xlab=xLabel, ylab=yLabel, 
                   col=colors, plot.grid=TRUE, grid.inc=9, horiz=TRUE, xpd=FALSE)

		if(i == 1) {		
			# The y-axis with labels for each group
			axis(2, labels=labels, at=mp, tck=0, cex.axis=0.75)
		}
		
		# The x-axis with the values going from xMin to xMax
		r <- range(valuesTmp, na.rm=TRUE)
		diffR <- r[2]-r[1]
		
		if(diffR > 25) {
			byVal <- 25
		} else if(diffR > 10 & diffR < 25) {
			byVal <- 5
		} else {
			byVal <- 0.5
		}
		
		axis(1, at=seq(xMin, xMax, by=byVal))	
	}
  
	if(verbose) {
		str(opar)		
	}
	
	par(opar)
}


