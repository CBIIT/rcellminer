#' Make a simple 2d plot using two variables with ggplot2 
#' 
#' @param df a data.frame with at least two columns 
#' @param xCol the name of the column in df with the "x" data
#' @param yCol the name of the column in df with the "y" data 
#' @param xPlotLabel the x plot label
#' @param yPlotLabel the y plot label
#' @param title the plot title, if null the correlation will appear (DEFAULT: null)
#' @param showColorTissues boolean, whether to show tissue colors 
#' @param showLegend boolean, whether to show the legend 
#' @param showTrendLine boolean, whether to show the trendline 
#' @param showTitle NNNNN
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
#' plotCellMiner2D(df, xCol="x", yCol="y", xPlotLabel="SLFN11", yPlotLabel="94600")
#' plotCellMiner2D(df, xCol="x", yCol="y", showTrendLine = FALSE, showTitle = FALSE)
#' plotCellMiner2D(df, xCol="x", yCol="y", showTrendLine = FALSE, showLegend = FALSE)
#' }
#' 
#' @author Augustin Luna <augustin AT mail.nih.gov>
#' 
#' @importFrom ggplot2 ggplot geom_point theme_bw scale_colour_manual xlab ylab geom_smooth ggtitle theme
#' 
#' @concept rcellminer
#' @export
plotCellMiner2D <- function(df, xCol="x", yCol="y", xPlotLabel=xCol, yPlotLabel=yCol, 
														title=NULL, showColorTissues=TRUE, showLegend=TRUE, 
														showTrendLine=TRUE, showTitle=TRUE) {
	if(is.null(title)) {
		corResults <- cor.test(df[,xCol], df[,yCol], use="pairwise.complete.obs")
		title <- paste0(paste(yPlotLabel, '~', xPlotLabel),
										', r=', round(corResults$estimate, 2),
										' p=', signif(corResults$p.value, 2))
	}

	# Get colors 
	t1 <- loadNciColorSet(returnDf = TRUE)
	t2 <- unique(t1[, c("tissues", "colors")])
	colors <- t2$colors
	names(colors) <- t2$tissues
	
	# Plot image
	p1 <- ggplot(data=df, aes_string(x=xCol, y=yCol))
	p1 <- p1 + geom_point(aes(colour = tissues))
	p1 <- p1 + theme_bw()
	
	if(showColorTissues) {
		p1 <- p1 + scale_colour_manual(values=colors)
	}

	if(!is.null(xPlotLabel)) {
		p1 <- p1 + xlab(xPlotLabel)		
	} else {
		p1 <- p1 + xlab(xCol)				
	}

	if(!is.null(yPlotLabel)) {
		p1 <- p1 + ylab(yPlotLabel)		
	} else {
		p1 <- p1 + ylab(yCol)				
	}
	
	if(showTrendLine) {
		p1 <- p1 + geom_smooth(method = "lm", se = FALSE, color = "red")
	}
	
	if(showTitle) {
		p1 <- p1 + ggtitle(title)
	}
	
	if(!showLegend) {
		p1 <- p1 + theme(legend.position="none")
	}
	
	return(p1)
}