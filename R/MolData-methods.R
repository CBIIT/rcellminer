#' Returns a MolData object.
#'
#' @param eSetList A list of eSet objects for a common set of samples.
#' @param sampleData A MIAxE object capturing sample and other data set information.
#' @param .Object An object: see "new()" documentation in "methods" package. 
#' 
#' @return A MolData object.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod initialize
setMethod(f="initialize",
					signature="MolData",
					definition=function(.Object,
															eSetList,
															sampleData){
						.Object@eSetList <- eSetList
						.Object@sampleData <- sampleData
						
						if (validObject(.Object)){
							return(.Object)
						}
					})


#' Returns a MolData object.
#'
#' @param eSetList A list of eSet objects for a common set of samples.
#' @param sampleData A MIAxE object capturing sample and other data set information.
#' @param ... Other possible parameters.
#' 
#' @return A MolData object.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("MolData",
					signature(eSetList="list", sampleData="MIAxE"),
					function(eSetList, sampleData, ...){
						.MolData(eSetList, sampleData)
					})

#' Returns a list of eSet objects.
#'
#' @param object MolData object for which a list of eSet objects is to be returned.
#' @return A list of eSet objects.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getESetList", signature = "MolData", function(object) { object@eSetList })


#' Returns an indexed eSet object from a MolData object eSet list.
#'
#' @param x A MolData object.
#' @param i Index or named item in MolData object eSet list.
#' @return An indexed eSet object from a MolData object eSet list.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("[[", signature = "MolData", function(x, i) { getESetList(x)[[i]] })


#' Assigns an eSet object to a specified position in a MolData object eSet list.
#'
#' @param x A MolData object.
#' @param i Index or named item in MolData object eSet list.
#' @param value An eSet object to be assigned.
#' @return An eSet object to a specified position in a MolData object eSet list.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("[[<-", signature(x="MolData"), function(x, i, value) {
	x@eSetList[[i]] <- value
	if (validObject(x)){
		return(x)
	}
})


#' Returns a list of feature data matrices.
#'
#' @param object MolData object for which a list of feature data matrices is to be returned.
#' @return A list of feature data matrices.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getAllFeatureData", signature = "MolData", function(object) {
	lapply(getESetList(object), exprs) })

#' Returns a data frame with sample information.
#'
#' @param object MolData object for which sample data is to be returned.
#' @return A data frame with sample information.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getSampleData", signature = "MolData", function(object) {
	as.data.frame(do.call(cbind, samples(object@sampleData)), stringsAsFactors=FALSE)
})

#' Returns a list of data frames with feature information.
#'
#' @param object MolData object for which feature data is to be returned.
#' @return A named list of data frames with feature information for 
#'  available molecular data types.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setMethod
setMethod("getFeatureAnnot", signature = "MolData", function(object) {
	featureDat <- lapply(getESetList(object), FUN = function(x) { 
			as(featureData(x), "data.frame") 
		})
	return(featureDat)
})