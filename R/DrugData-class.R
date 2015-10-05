#' An S4 class to represent drug activity and related data recorded for a set of
#' biological samples.
#'
#' @slot act An eSet object containing drug activity data across a set of biological samples.
#' @slot repeatAct An eSet object containing repeat drug activity experiment data with
#' respect to the same samples associated with act.
#' @slot sampleData A MIAxE object capturing sample and other data set information.
#' 
#' @param ... Other possible parameters.
#'
#' @exportClass DrugData
#' @concept rcellminerData
#' 
#' @importFrom methods setClass
#' @importFrom Biobase samples
.DrugData <- setClass(Class="DrugData",
										 slots=c(act = "eSet",
										 				repeatAct = "eSet",
										 				sampleData = "MIAxE"),
										 validity=function(object){
										 	sampleNames <- samples(object@sampleData)[[1]]
										 	
										 	if (!identical(colnames(exprs(object@act)), sampleNames)){
										 		return("eSet act samples must match samples in sampleData.")
										 	}
										 	
										 	if (!identical(colnames(exprs(object@repeatAct)), sampleNames)){
										 		return("eSet repeatAct samples must match samples in sampleData.")
										 	}
										 	
										 	return(TRUE)
										 })



