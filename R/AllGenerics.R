#--------------------------------------------------------------------------------------------------
# GENERIC FUNCTION DEFINITIONS
#--------------------------------------------------------------------------------------------------

#' Returns a DrugData object.
#'
#' @param act An eSet object containing drug activity data across a set of biological samples.
#' @param repeatAct An eSet object containing repeat drug activity experiment data with
#' respect to the same samples associated with act.
#' @param sampleData A MIAxE object capturing sample and other data set information.
#' @param ... Other possible parameters.
#' @return A DrugData object.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric("DrugData", 
					 function(act, repeatAct, sampleData, ...) { standardGeneric("DrugData") } )


#' Returns a MolData object.
#'
#' @param eSetList A list of eSet objects for a common set of samples.
#' @param sampleData A MIAxE object capturing sample and other data set information.
#' @param ... Other possible parameters.
#' @return A MolData object.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric("MolData", 
					 function(eSetList, sampleData, ...) { standardGeneric("MolData") } )


#' Returns a list of feature data matrices.
#'
#' @param object Object for which a list of feature data matrices is to be returned.
#' @param ... Other possible parameters.
#' @return A list of feature data matrices.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric(name="getAllFeatureData",
           def=function(object, ...) { standardGeneric("getAllFeatureData") })

#' Returns a list of eSet objects.
#'
#' @param object Object for which a list of eSets is to be returned.
#' @param ... Other possible parameters.
#' @return A list of eSet objects.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric(name="getESetList",
					 def=function(object, ...) { standardGeneric("getESetList") })

#' Returns a data frame with sample information.
#'
#' @param object Object for which sample data is to be returned.
#' @param ... Other possible parameters.
#' @return A data frame with sample information.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric(name="getSampleData",
           def=function(object, ...) { standardGeneric("getSampleData") })


#' Returns a list of data frames with feature information.
#'
#' @param object Object for which feature data is to be returned.
#' @param ... Other possible parameters.
#' @return A list of data frames with feature information.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric(name="getFeatureAnnot",
					 def=function(object, ...) { standardGeneric("getFeatureAnnot") })


#' Returns an eSet object with drug activity data.
#'
#' @param object Object for which drug activity data is to be returned.
#' @param ... Other possible parameters.
#' @return An eSet object with drug activity data.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric(name="getAct",
           def=function(object, ...) { standardGeneric("getAct") })


#' Returns an eSet object with drug repeat activity experiment data.
#'
#' @param object Object for which drug repeat activity experiment data is to be returned.
#' @param ... Other possible parameters.
#' @return An eSet object with drug repeat activity experiment data.
#'
#' @concept rcellminerData
#' @export
#' 
#' @importFrom methods setGeneric
setGeneric(name="getRepeatAct",
           def=function(object, ...) { standardGeneric("getRepeatAct") })

#--------------------------------------------------------------------------------------------------
