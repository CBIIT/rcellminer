#' Make a simple 2d plot using two variables with ggplot2 
#' 
#' @param df a data.frame with at least two columns 
#' @param xCol the name of the column in df with the "x" data. See Note
#' @param yCol the name of the column in df with the "y" data. See Note
#' @param xLabel the x plot label
#' @param yLabel the y plot label
#' @param title the plot title, if null the correlation will appear (DEFAULT: NULL)
#' @param classCol the name of the column with the classes. Values in column of df must be a factor (DEFAULT: NULL)
#' @param tooltipCol the name of the column used for tooltips when plotted with plotly
#' @param colorPalette a named vector with the names classes and value colors (DEFAULT: NULL)
#' @param showLegend boolean, whether to show the legend (DEFAULT: FALSE)
#' @param showTrendLine boolean, whether to show the trendline 
#' @param showTitle boolean, whether to show the title
#' @param singleColor a color to be used for all points when a color palette is not provided (DEFAULT: blue)
#' @param alpha value from 0-1, where 0 indicates transparent points (DEFAULT: 1, not transparent)
#' @param numberColPrefix a prefix to add to column names that start with a number that causes issues with ggplot (DEFAULT: X)
#' @param xLimVal a two entry vector (min, max) to set the x-axis
#' @param yLimVal a two entry vector (min, max) to set the y-axis
#' @param pointSize size of points on plot (DEFAULT: 3)
#' 
#' @note TROUBLESHOOTING NOTES: 1) Avoid ":" in colnames
#' 
#' @return a ggplot object
#' 
#' @examples
#' \dontrun{
#' # Load data
#' nci60DrugActZ <- exprs(getAct(rcellminerData::drugData))
#' nci60GeneExpZ <- getAllFeatureData(rcellminerData::molData)[["exp"]]
#' # Load colors
#' colorTab <- loadNciColorSet(returnDf=TRUE)
#' tissueColorTab <- unique(colorTab[, c("tissues", "colors")])
#' # Merge data
#' df <- as.data.frame(t(rbind(nci60DrugActZ["94600",], nci60GeneExpZ["SLFN11",])))
#' colnames(df) <- c("y", "x")
#' df <- cbind(df, colorTab)
#' # Plot data
#' plotCellMiner2D(df, xCol="x", yCol="y", xLabel="SLFN11", yLabel="94600")
#' plotCellMiner2D(df, xCol="x", yCol="y", showTrendLine = FALSE, showTitle = FALSE)
#' plotCellMiner2D(df, xCol="x", yCol="y", showTrendLine = FALSE, showLegend = FALSE)
#' }
#' 
#' @author Augustin Luna <augustin AT mail.nih.gov>
#' 
#' @note Uses ggplot aes_string() which uses parse() to turn your text expression into a proper R symbol that can be resolved within the data.frame. Avoid numbers and spaces in 
#' 
#' @importFrom ggplot2 ggplot geom_point theme_bw scale_colour_manual xlab ylab geom_smooth ggtitle theme aes_string
#' 
#' @concept rcellminer
#' @export
plotCellMiner2D <- function(df, xCol="x", yCol="y", xLabel=xCol, yLabel=yCol, 
														title=NULL, colorPalette=NULL, classCol=NULL, tooltipCol=NULL, 
														showLegend=FALSE, showTrendLine=TRUE, showTitle=TRUE, singleColor="#0000FF",
														alpha=1, numberColPrefix="X", xLimVal=NULL, yLimVal=NULL, pointSize=3) {
	
	# nci60DrugActZ <- exprs(getAct(rcellminerData::drugData))
	# nci60GeneExpZ <- getAllFeatureData(rcellminerData::molData)[["exp"]]
	# # Load colors
	# colorTab <- loadNciColorSet(returnDf=TRUE)
	# tissueColorTab <- unique(colorTab[, c("tissues", "colors")])
	# # Merge data
	# xCol <- "SLFN11"
	# yCol <- "94600"
	# classCol <- "tissues"
	# xLabel <- xCol
	# yLabel <- yCol
	# title <- NULL
	# df <- data.frame(y=nci60DrugActZ[yCol,], x=nci60GeneExpZ[xCol,])
	# yCol <- "X94600" # MUST NOT BE A NUMBER
	# colnames(df) <- c(yCol, xCol)
	# df <- cbind(df, colorTab)
	# df[, classCol] <- as.factor(df[, classCol])
	# colorPalette <- tissueColorTab[, "colors"]
	# names(colorPalette) <- tissueColorTab[, classCol]
	# colors <- rep("blue", nrow(df))
	# colors[1:10, "colors"] <- "red"
	# showLegend <- FALSE
	# showTrendLine <- TRUE
	# showTitle <- TRUE
	# alpha <- 1
	
	# Fix column names if they start with a number 
	if (grepl("^[0-9]", xCol)) {
		tmpXCol <- paste0(numberColPrefix, xCol)
		df[, tmpXCol] <- df[, xCol]
		xCol <- tmpXCol
	}

	if (grepl("^[0-9]", yCol)) {
		tmpYCol <- paste0(numberColPrefix, yCol)
		df[, tmpYCol] <- df[, yCol]
		yCol <- tmpYCol
	}

	# Create title 
	if(is.null(title)) {
		corResults <- cor.test(df[,xCol], df[,yCol], use="pairwise.complete.obs")
		title <- paste0(paste(yLabel, 'vs.', xLabel),
										'\nPearson correlation (r)=', round(corResults$estimate, 2),
										', p-value=', formatC(signif(corResults$p.value, 2)))
	}

	# Plot image
	p1 <- ggplot(data=df, aes_string(x=xCol, y=yCol))
	p1 <- p1 + theme_bw()
	
	if(!is.null(colorPalette) && !is.null(classCol)) {
		p1 <- p1 + suppressWarnings(geom_point(aes_string(color=classCol, text=tooltipCol), alpha=alpha, size = pointSize))
		p1 <- p1 + scale_colour_manual(name="", values=colorPalette)
	} else {
		p1 <- p1 + suppressWarnings(geom_point(aes_string(text=tooltipCol), color=singleColor, alpha=alpha, size = pointSize))
	}

	if(!is.null(xLabel)) {
		p1 <- p1 + xlab(xLabel)		
	} else {
		p1 <- p1 + xlab(xCol)				
	}

	if(!is.null(yLabel)) {
		p1 <- p1 + ylab(yLabel)		
	} else {
		p1 <- p1 + ylab(yCol)				
	}
	
	if(showTrendLine) {
		p1 <- p1 + geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.25)
	}
	
	if(showTitle) {
		p1 <- p1 + ggtitle(title)
	}
	
	if(!showLegend) {
		p1 <- p1 + theme(legend.position="none", plot.title = element_text(size=14),
										 axis.title.x = element_text(size=20),
										 axis.title.y = element_text(size=20))
	}
	
	if(!is.null(xLimVal)) {
		p1 <- p1 + xlim(xLimVal[1], xLimVal[2])
	}
	
	if(!is.null(yLimVal)) {
		p1 <- p1 + ylim(yLimVal[1], yLimVal[2])
	}
	
	return(p1)
}


