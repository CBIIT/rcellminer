#' Select features that are correlated with a given feature (or one or more features
#' from a set of features). 
#' 
#' @param Y a vector or matrix; rows from X will be correlated with Y if Y is a vector
#' or with rows of Y, if Y is a matrix.
#' @param X a matrix of values that will be compared with Y (vector) or rows of Y (matrix)
#' @param corThreshold the minimum correlation threshold for the row to be returned
#' @param useAbsCor a logical value indicating whether absolute correlations should be used
#' (default=TRUE).
#' @return a matrix of rows of X correlated with Y (if Y is a vector) or correlated with
#' at least one row of Y if Y is a matrix or data frame. 
#' 
#' @examples 
#' vec <- runif(10)
#' mat <- matrix(runif(100), 10, 10)
#' selectCorrelatedRows(vec, mat)
#' 
#' @concept rcellminer
#' @export
selectCorrelatedRows <- function(Y, X, corThreshold=0.10, useAbsCor=TRUE){
	crossCorMat <- crossCors(Y, X)$cor
	if (useAbsCor){
		crossCorMat <- abs(crossCorMat)
	}
	
	if (is.vector(Y) || (nrow(Y) ==1)){
		#decCorOrder <- order(crossCorMat, decreasing=TRUE)
		#X <- X[decCorOrder, ]
		#crossCorMat <- crossCorMat[, decCorOrder, drop=FALSE]
		#stopifnot(identical(rownames(X), colnames(crossCorMat)))
		return(X[which(crossCorMat > corThreshold), , drop=FALSE])
	}
	else {
		return(X[apply(X=(crossCorMat > corThreshold), MARGIN=2, FUN=any), , drop=FALSE])
	}
}
