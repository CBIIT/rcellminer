# NOTE: Size is not automatically set for rChartsAlternative output
plotHeight <- 900
plotWidth <- 900

tooltipCol <- "tooltip"

isDrugActivityDataType <- function(prefix){
	# TO DO: Make configurable.
	drugActTypePrefixes <- "act"
	if (prefix %in% drugActTypePrefixes){
		return(TRUE)
	} else{
		return(FALSE)
	}
}

isGeneProtDataType <- function(prefix){
	# TO DO: Make configurable.
	geneProtDataTypePrefixes <- c("cop", "mut", "met", "exp", "xai", "swa", "pro")
	if (prefix %in% geneProtDataTypePrefixes){
		return(TRUE)
	} else{
		return(FALSE)
	}
}