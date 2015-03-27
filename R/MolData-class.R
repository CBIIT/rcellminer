#' An S4 class to represent molecular data recorded for a set of biological samples.
#'
#' @slot eSetList A list of eSet objects for a common set of samples.
#' @slot sampleData A MIAxE object capturing sample and other data set information.
#' 
#' @param ... Other possible parameters.
#'
#' @exportClass MolData
#' @concept rcellminerData
#' 
#' @importFrom methods setClass
.MolData <- setClass(Class="MolData",
										slots=c(eSetList = "list",
														sampleData = "MIAxE"),
										validity = function(object){
											if (length(getESetList(object)) > 0){
												
												if (!all(vapply(getESetList(object), function(x) is(x, "eSet"), logical(1)))){
													return("All added elements must be subclasses of eSet.")
												}
												
												sampleNames <- samples(object@sampleData)[[1]]
												eSetSamplesMatchSampleData <- vapply(getESetList(object),
																														 function(x) identical(colnames(exprs(x)), 
																														 											sampleNames), logical(1))
												
												if (!all(eSetSamplesMatchSampleData)){
													return("All eSet samples must match samples in sampleData.")
												}
												
											}
											
											return(TRUE)
										})

