#' Description: Produces CellMiner-like plots in R 
#' 
#' @param drugAct a matrix of drug activity values (cell lines as columns, drug entries as rows)
#' @param molData a list of matricies a molecular 
#' @param plots a vector of characters denoting the plots to include and the order (e.g. c("mut", "drug", "cop"). 
#'   Currently, supported entries mutations (mut), drug activities (drug), copy number variations (cop)
#' @param nsc a string NSC ID that will be plotted when a "drug" entry appears in the plots vector  
#' @param gene a string HUGO gene symbol for which the "mut", "cop", or "exp" plots will be produced if in plots vector 
#' @param sub a vector of strings with sub-titles for each plot
#' @param xLimits a 2 number vector with the the minimum and maximum X-axis values (default: c(-3, 3)) 
#' @param xLabel a string for the default X-axis label
#' @param verbose a boolean to show debugging information
#' @return None
#' 
#' @author Augustin Luna <augustin AT mail.nih.gov>
#' 
#' @importFrom gplots barplot2
#' 
#' @examples
#' drugAct <- exprs(getAct(rcellminerData::drugData))
#' molDataMats <- getMolDataMatrices()
#' plotDataMats <- molDataMats[c("exp", "cop", "mut")]
#' plotCellMiner(drugAct, plotDataMats, 
#' plots=c("mut", "drug", "cop"), nsc="94600", gene="TP53")
#' 
#' @concept rcellminer
#' @export
plotCellMiner <- function(drugAct, molData, plots, nsc=NULL, gene=NULL, sub=NULL, 
                          xLimits=c(-3, 3), xLabel="Z-scores", verbose=FALSE) {
  # Make sure the nsc is a character string otherwise NSCs may be treated as numeric indicies
  nsc <- as.character(nsc)
  
  expPrefix <- "exp"
  mutPrefix <- "mut"
  copPrefix <- "cop"
  
	expData <- molData[[expPrefix]]
	mutData <- molData[[mutPrefix]]
	copData <- molData[[copPrefix]]
	
  drugIdx <- 1
  expIdx <- 1
	mutIdx <- 1
	copIdx <- 1

  mainTitle <- NULL

  if(!is.null(sub)) {
    subTitle <- sub
  } else {
    subTitle <- rep(NULL, length(plots))
  }
  
	values <- NULL

	for(p in plots) {
		if(p == "drug") {
			#mainTitle <- c(mainTitle, paste("NSC # ", nsc[drugIdx], " average Z scores", collapse=""))
			mainTitle <- c(mainTitle, paste("NSC", nsc[drugIdx], collapse=""))
			values <- cbind(values, as.numeric(drugAct[nsc[drugIdx],]))
      drugIdx <- drugIdx + 1
		}
		
		if(p == expPrefix) {
			mainTitle <- c(mainTitle, paste(gene[expIdx], " expression", collapse=""))
			values <- cbind(values, as.numeric(expData[paste0(expPrefix, gene[expIdx]),]))
			expIdx <- expIdx + 1
		}

		if(p == mutPrefix) {
			mainTitle <- c(mainTitle, paste(gene[mutIdx], " mutation", collapse=""))
			values <- cbind(values, as.numeric(mutData[paste0(mutPrefix, gene[mutIdx]),]))
			mutIdx <- mutIdx + 1
		}

		if(p == copPrefix) {
			mainTitle <- c(mainTitle, paste(gene[copIdx], " copy number", collapse=""))
			values <- cbind(values, as.numeric(copData[paste0(copPrefix, gene[copIdx]),]))
			copIdx <- copIdx + 1
		}	
	}

	labels <- colnames(drugAct) 

	### Initiate Graphic  
	# NOTE: To make longer gene names fit use cex.axis 
	# in Heatmap section to scale down font or MAYBE increase the 
	# left margin c(BOTTOM, LEFT, TOP, RIGHT) in Heatmap and GI50
	
  # Save default for resetting
  opar <- par(no.readonly = TRUE)
  
	figures <- length(plots) 
	figureCols <- 3
	m <- NULL

	for(i in 1:figures) {
		tmp <- NULL
		
		for(j in 1:figureCols) {
			tmp <- cbind(tmp, rep(i, 10))
		}

		m <- cbind(m, tmp) 
	}
	layout(m, heights = lcm(12)) #widths = lcm(4), 
	#layout.show() 

	# Labels
	xLabel <- xLabel
	yLabel <- ""

	# Limits 
	xMin <- xLimits[1]
	xMax <- xLimits[2]
	
	par(las=1, cex=0.5)
	# Margins: Bottom, Left, Top, Right
	par(oma=c(4,4,4,4)) 
	
	for(i in 1:figures) {
		if(verbose) {
			str(figures) 
			str(values) 
			str(labels) 
			str(values[,i]) 
			str(loadNciColorSet())
		}
		
		# Match CellMiner plots with BR lines at the top
	  labels <- rev(labels) 
	  valuesTmp <- rev(values[,i])
	  colors <- rev(loadNciColorSet())
	  
	  #plotNci60Profile(values[,i], xMin, xMax, mainTitle[i], xLabel, yLabel, as.character(colors), labels)  
                     
		# Return the midpoints 
		mp <- barplot2(valuesTmp, axes=FALSE, axisnames=FALSE, xlim=c(xMin, xMax), 
                   main=mainTitle[i], sub=subTitle[i], xlab=xLabel, ylab=yLabel, 
                   col=as.character(colors), plot.grid=TRUE, grid.inc=9, horiz=TRUE, xpd=FALSE)

		if(i == 1) {		
			# The y-axis with labels for each group
			axis(2, labels=labels, at=mp, tck=0, cex.axis=0.75)
		}
		
		# The x-axis with the values going from xMin to xMax
		axis(1, at=seq(xMin, xMax, by=1))	
	}
  
	if(verbose) {
		str(opar)		
	}
	
	par(opar)
}


