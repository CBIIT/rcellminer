#' Extract from a list of matrices the data associated with a set of features.
#' 
#' @param featureSet a character vector of feature names.
#' @param dataMatList a list of matrices with feature data organized along the rows,
#' and feature names accessible via rownames(dataMatList).
#' @param excludeMissingFeatures a logical value indicating whether features whose data 
#' cannot be found in any matrices in dataMatList should be excluded in the output.
#' (default=TRUE).
#' @return a single matrix containing data for all features in featureSet.
#' 
#' @examples 
#' featureSet <- c("expSLFN11", "mutSLX4")
#' molDataMats <- getMolDataMatrices()
#' featureData <- getFeatureDataFromMatList(featureSet, molDataMats)
#' 
#' @concept rcellminer
#' @export
getFeatureDataFromMatList <- function(featureSet, dataMatList, excludeMissingFeatures=TRUE){
  iMax <- which.max(vapply(dataMatList, ncol, integer(1)))
  numCols <- ncol(dataMatList[[iMax]])
  colNames <- colnames(dataMatList[[iMax]])
  featureSetMat <- matrix(NA, nrow=length(featureSet), ncol=numCols)
  rownames(featureSetMat) <- featureSet
  colnames(featureSetMat) <- colNames
  
  for (feature in featureSet){
    iMatchedMat <- which(vapply(dataMatList, 
                                function(mat) is.element(feature, rownames(mat)),
                                logical(1)))
    if (length(iMatchedMat) == 0){
      next
    }
    iMatchedMat <- iMatchedMat[1]
    featureSetMat[feature, ] <- dataMatList[[iMatchedMat]][feature, , drop=FALSE]
  }
  
  if (excludeMissingFeatures){
    rowHasData <- apply(X=featureSetMat, MARGIN=1, FUN=function(x) !all(is.na(x)))
    featureSetMat <- featureSetMat[rowHasData, , drop=FALSE]
  }
  
  return(featureSetMat)
}