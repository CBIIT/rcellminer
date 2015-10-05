#' Produces a barplot of the average values for a set of NSCs with a error bar (one standard deviation)
#' 
#' @param drugAct a matrix of drug activity values (cell lines as columns, drug entries as rows)
#' @param drugs a vector of NSC IDs whose values will be averaged by cell line
#' @param mainLabel a main label for the plot
#' @param pdfFilename a string file name for a PDF plot, no file output will be produced if this is not provided 
#' @param statistic a string, either 'mean' or 'median' (Default: mean)
#' 
#' @return no values are returned
#' 
#' @examples 
#' drugAct <- exprs(getAct(rcellminerData::drugData))
#' drugs <- rownames(drugAct)[1:8]
#' plotDrugSets(drugAct, drugs, "Test")
#' 
#' @importFrom gplots barplot2
#' 
#' @concept rcellminer
#' @export
plotDrugSets <- function(drugAct, drugs, mainLabel="", pdfFilename=NULL, statistic="mean") {	
	if (!is.data.frame(drugAct)){
		drugAct <- as.data.frame(drugAct)
	}
	data <- drugAct[drugs,]

	if(statistic == "mean") {
		tmpData <- sapply(data, mean, na.rm=TRUE)
	} else if (statistic == "median") {
		tmpData <- sapply(data, median, na.rm=TRUE)
	} else {
		stop("Unknown statistic. Please choose 'mean' or 'median'")
	}

	stDevs <- sapply(data, sd, na.rm=TRUE)

	labels <- colnames(drugAct)
	#mainLabel <- paste(annotation, "Total drugs:", length(drugs), sep=" ")

	# Limits 
	yMax <- 4
	yMin <- -4

	if(!is.null(pdfFilename)) {
		pdf(pdfFilename)
	}

	#NCI60 color set
	#colors <- c(rep("FF0000", 8), rep("FFFF00", 2), rep("610B0B", 7), rep("0489B1", 9), rep("4B610B", 9), rep("86B404", 7), rep("DF7401", 7), rep("8A4B08", 6), rep("0B4C5F", 5))
	#colors <- c(rep("red", 8), rep("yellow", 2), rep("brown4", 7), rep("steelblue3", 9), rep("chartreuse4", 9), rep("olivedrab3", 7), rep("orange1", 7), rep("orange4", 6), rep("dodgerblue4", 5))
	#colors <- rev(colors)
	colors <- loadNciColorSet()
	
	# Return the midpoints 
	mp <- barplot2(tmpData, axes=FALSE, axisnames=FALSE, ylim=c(yMin, yMax), main=mainLabel, ylab="Z-Scores", col=as.character(colors), plot.grid=TRUE, grid.inc=9, horiz=FALSE)

	par(las=3, cex=0.5)
	# The x-axis with labels for each group
	axis(1, labels=labels, at=mp)
	# The y-axis with the values going from yMin to yMax
	axis(2, at=seq(yMin, yMax, by=1))

	# Get standard deviation of each group
	stDevsMat <- matrix(stDevs, 2)
	# Plot the vertical lines of the error bars
	segments(mp, tmpData - stDevs, mp, tmpData + stDevs, lwd=2)
	# Now plot the horizontal bounds for the error bars
	#segments(mp - 0.1, tmpData - stDevs, mp + 0.1, tmpData - stDevs)
	#segments(mp - 0.1, tmpData + stDevs, mp + 0.1, tmpData + stDevs)
	
	if(!is.null(pdfFilename)) {
		dev.off()
	}
}
