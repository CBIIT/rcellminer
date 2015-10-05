#' Returns data to plot CellMiner plots
#' 
#' @param returnDf a boolean if a data.frame with all information (default: FALSE)
#' @return a vector of colors as strings or a data.frame with dataType, label, xMin, xMax
#'   
#' @examples
#' loadCellminerPlotInfo()  
#'   
#' @concept rcellminer
#' @export
loadCellminerPlotInfo <- function(returnDf=FALSE) {
	featureTypes <- c("drug", "cop", "exp", "xai", "exo", "mut", "mir", "pro", "mda", 
								   "mdaage", "mdais_epithelial", "mdais_p53_mut", "mdamdr", "mdadoublingtime")
	
	labels <- c("Drug Activity (Z-Score)",
						 "Copy Number Deviation\n(0 is Copy Number 2)", 
						 "Transcript Intensity (Z-Score)", 
						 "Transcript Intensity (log2)",
						 "% Variant Allele", 
						 "Deleterious Mutations Present (Binary)", 
						 "miRNA Abundance (log2)",
						 "Protein Abundance (log2)",
						 NA,
						 "Age (Years)",
						 "Is Epithelial (Binary)",
						 "Is TP53 Mutated (Binary)",
						 "Assay Units (Flow Cytometry)",
						 "Doubling Time (Hours)")
	
	titlePrefixes <- c("NSC: ", "Copy Number: ", "Expression: ", "Expression: ", 
										 "Allele: ", "Mutations: ", "miRNA: ", "Protein: ", "Metadata: ", 
										 "Age: ", "Epithelial: ", "TP53 Mutations: ", "MDR: ", "Doubling Time: ")
	
	xMin <- c(-3, NA, -3, NA, 0, 0, NA, NA, NA, NA, 0, 0, NA, NA)
	xMax <- c(3, NA, 3, NA, 100, 1, NA, NA, NA, NA, 1, 1, NA, NA)
	
	results <- data.frame(featureType=featureTypes, label=labels, xMin=xMin, xMax=xMax, 
												titlePrefix=titlePrefixes, stringsAsFactors=FALSE)
	
	if(returnDf) {
		return(results)
	} else {
		return(labels)     
	}
}
