#' Returns a DrugData object.
#' 
#' @param act An eSet object containing drug activity data across a set of biological samples.
#' @param repeatAct An eSet object containing repeat drug activity experiment data with
#' respect to the same samples associated with act.
#' @param sampleData A MIAxE object capturing sample and other data set information.
#' @param .Object An object: see "new()" documentation in "methods" package. 
#' 
#' @return A DrugData object.
#'
#' @concept rcellminerData
#' @export
#' @note Seems to be required for definition of a constructor.
#' 
#' @importFrom methods setMethod initialize
setMethod(f="initialize",
					signature="DrugData",
					definition=function(.Object, act, repeatAct, sampleData){
						.Object@act <- act
						.Object@repeatAct <- repeatAct
						.Object@sampleData <- sampleData
						
						if (validObject(.Object)){
							return(.Object)
						}
					})



#' Returns a DrugData object.
#'
#' @param act An eSet object containing drug activity data across a set of biological samples.
#' @param repeatAct An eSet object containing repeat drug activity experiment data with
#' respect to the same samples associated with act.
#' @param sampleData A MIAxE object capturing sample and other data set information.
#' @param ... Other possible parameters.
#' 
#' @return A DrugData object.
#' 
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("DrugData",
					signature(act="eSet", repeatAct="eSet", sampleData="MIAxE"),
					function(act, repeatAct, sampleData, ...){
						.DrugData(act, repeatAct, sampleData)
					})


#' Returns a data frame with sample information.
#'
#' @param object DrugData object for which sample data is to be returned.
#' @return A data frame with sample information.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getSampleData", signature = "DrugData", function(object) {
	as.data.frame(do.call(cbind, samples(object@sampleData)), stringsAsFactors=FALSE)
})

#' Returns a list of data frames with feature information.
#'
#' @param object DrugData object for which feature data is to be returned.
#' @return A named list of data frames with feature information for drugs
#'  and drug repeat experiments.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getFeatureAnnot", signature = "DrugData", function(object) {
	featureDat <- list()
	featureDat$drug <- as(featureData(getAct(object)), "data.frame")
	featureDat$drugRepeat <- as(featureData(getRepeatAct(object)), "data.frame")
	return(featureDat)
})

#' Returns an eSet object with drug activity data.
#'
#' @param object DrugData object for which drug activity data is to be returned.
#' @return An eSet object with drug activity data.
#' 
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getAct", signature = "DrugData", function(object) {
	object@act
})

#' Returns an eSet object with drug repeat activity experiment data.
#'
#' @param object DrugData object for which drug repeat activity experiment data
#' is to be returned.
#' @return An eSet object with drug repeat activity experiment data.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getRepeatAct", signature = "DrugData", function(object) {
	object@repeatAct
})



