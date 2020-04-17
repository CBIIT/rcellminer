#' Computes the relative standard deviation values with respect to the columns of a matrix or data.frame.
#' 
#' @param dat a matrix or data.frame with numeric values.
#' @param onlyReturnMedian a logical value indicating whether only the median column RSD value should be returned (vs. all RSD values).
#' 
#' @return median RSD value over the data set columns or all RSD values, depending on value of onlyReturnMedian (default=TRUE).
#' 
#' @examples 
#' A <- matrix(rnorm(10*60), nrow=10)
#' getRsd(A)
#' getRsd(A, onlyReturnMedian=FALSE)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stats sd
getRsd <- function(dat, onlyReturnMedian=TRUE){
  rsdVals <- abs(apply(dat, MARGIN=2, FUN=sd, na.rm=TRUE) / apply(dat, MARGIN=2, FUN=mean, na.rm=TRUE))
  
  if (onlyReturnMedian){
    return(median(rsdVals, na.rm=TRUE))
  } else{
    return(rsdVals)
  }
}