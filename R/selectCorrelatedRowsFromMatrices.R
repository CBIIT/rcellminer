#' Select features that are correlated with a given feature (or one or more features
#' from a set of features), merging results from multiple candidate feature matrices.
#' 
#' @param Y a vector or matrix; rows from each matrix element of X will be correlated with 
#' Y if Y is a vector or with rows of Y, if Y is a matrix.
#' @param XList a list of matrices whose rows will be correlated with Y (vector) or rows of Y (matrix).
#' The rownames in each matrix element of XList must be specified to values that are unique with
#' respect to the total set of rownames (as derived from all matrices in XList).
#' @param corThreshold the minimum correlation threshold for the row to be returned
#' @param useAbsCor a logical value indicating whether absolute correlations should be used
#' (default=TRUE).
#' @return a matrix formed from rows of matrices in XList that are correlated with Y 
#' (if Y is a vector) or correlated with at least one row of Y if Y is a matrix or data frame.
#' 
#' @examples 
#' vec <- runif(10)
#' names(vec) <- 1:10
#' matList <- list(X1 = matrix(runif(100), 10, 10), X2 = matrix(runif(100), 10, 10))
#' rownames(matList$X1) <- paste0("X1_row_", 1:10)
#' colnames(matList$X1) <- paste0("X1_col_", 1:10)
#' rownames(matList$X2) <- paste0("X2_row_", 1:10)
#' colnames(matList$X2) <- paste0("X2_col_", 1:10)
#' selectCorrelatedRowsFromMatrices(vec, matList)
#' 
#' @concept rcellminer
#' @export
selectCorrelatedRowsFromMatrices <- function(Y, XList, corThreshold=0.10, useAbsCor=TRUE){  
	if (!is.list(XList)){
		stop("Xlist parameter must be set to a list of numeric matrices.")
	}
	
	# Run selectCorrelatedRows on Y and each matrix in XList; put results in corRowData matrix list.
	corRowData <- lapply(XList, FUN=function(x) selectCorrelatedRows(Y, x, corThreshold, useAbsCor))
	
	# Remove corRowData entries for matrices in XList that did not match Y (at corThreshold).
	corRowData <- corRowData[which(vapply(corRowData, FUN=nrow, FUN.VALUE=integer(1)) > 0)]
	
	# Gather all correlated feature names.
	featureNames <- unname(c(lapply(corRowData, FUN=rownames), recursive=TRUE))
	
	if (is.null(dim(Y))){
		# Y is a vector.
		numCols <- length(Y)
		colNames <- names(Y)
	} else{
		# Y is a matrix or data frame.
		numCols <- ncol(Y)
		colNames <- colnames(Y)
	}
	
	featureMat <- matrix(0, nrow=length(featureNames), ncol=numCols)
	colnames(featureMat) <- colNames
	
	if (nrow(featureMat) > 0){
		rownames(featureMat) <- featureNames
		
		for (matName in names(corRowData)){
			featureMat[rownames(corRowData[[matName]]), ] <- corRowData[[matName]] 
		}
	}
	
	return(featureMat)
}