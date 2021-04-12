#' Returns a vector of log activity range values for set of compounds.
#' 
#' @param nscSet a character vector specifying NSC identifier(s) for compound(s) of interest.
#' @param computeIQR logical value indicated whether inter-quartile range is to be computed;
#' otherwise absolute range is computed (default=FALSE).
#' @return a numeric vector of NCI-60 log activity (-logGI50) range values indexed 
#' by the identifiers in nscSet.
#' 
#' @examples 
#' nscSet <- c("609699", "740")
#' getDrugActivityRange(nscSet)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stats IQR
getDrugActivityRange <- function(nscSet, computeIQR=FALSE){
	drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
  
  nscSet <- as.character(nscSet)
  actRangeVals <- numeric(length(nscSet))
  names(actRangeVals) <- nscSet
  
  for (nscStr in nscSet){
    iNsc <- which(as.character(drugRepeatAnnot$NSC) == nscStr)
    if (length(iNsc) == 0){
      warning(paste("No -logGI50 activity available for compound ", nscStr, ".", sep = ""))
      actRangeVals[nscStr] <- NA_real_
      next
    }
    
    avgNegLogGI50Data <- getDrugActivityData(nscStr)
    
    if (computeIQR){
      actRangeVals[nscStr] <- IQR(as.numeric(avgNegLogGI50Data[nscStr, ]), na.rm = TRUE)
    } else{
      actRangeVals[nscStr] <- max(as.numeric(avgNegLogGI50Data[nscStr, ]), na.rm = TRUE) - 
      	min(as.numeric(avgNegLogGI50Data[nscStr, ]), na.rm = TRUE) 
    }
  }
  
  return(actRangeVals)
}