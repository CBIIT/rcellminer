#' Calculate quantile for the columns in a matrix
#' 
#' @param X the matrix
#' @param prob a numeric probablity 
#' @param naRm a boolean, whether to remove NAs
#' @param onlyNonzeroVals a boolean, whether to only include non-zero values
#' @return a vector of quantiles
#' 
#' @examples
#' getColumnQuantiles(matrix(1:25, nrow=5), prob = 0.5)
#' 
#' @concept rcellminer
#' @export
getColumnQuantiles <- function(X, prob, naRm=FALSE, onlyNonzeroVals=FALSE){
	colQuantiles <- vector(mode="numeric", length=ncol(X))
	for (j in seq(colQuantiles)){
		if (onlyNonzeroVals){
			selector <- X[,j] != 0
		} else {
			selector <- seq(X[,j])
		}
		colQuantiles[j] <- quantile(x=X[selector,j], probs=prob, na.rm=naRm)
	}
	return(colQuantiles)
}
