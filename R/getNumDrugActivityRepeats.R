#' Returns a vector indicating the number of drug activity repeat experiments with available data
#' for each member of a set of compounds.
#' 
#' @param nscSet a character vector specifying NSC identifier(s) for compound(s) of interest.
#' @param onlyCellMinerExps a logical value indicating whether to return only the number of
#' experiments with data included in CellMiner (default=TRUE).
#' @return a numeric vector, indexed by nscSet, indicating the number of drug activity 
#' repeat experiments for each one of its compounds. 
#' 
#' @examples 
#' nscSet <- c("1", "17", "89", "609699")
#' getNumDrugActivityRepeats(nscSet)
#'
#' @concept rcellminer
#' @export
getNumDrugActivityRepeats <- function(nscSet, onlyCellMinerExps=TRUE){
	drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
  
  nscSet <- as.character(nscSet)
  numExps <- integer(length(nscSet))
  names(numExps) <- nscSet
  
  for (nscStr in nscSet){
    iNsc <- which(drugRepeatAnnot$nsc == nscStr)
    if (length(iNsc) == 0){
      warning(paste("No replicate experiment data available for compound ", nscStr, ".", sep = ""))
      numExps[nscStr] <- NA_real_
      next
    }
    
    if (onlyCellMinerExps){
      numExps[nscStr] <- sum(drugRepeatAnnot$used_in_zscore[iNsc])
    } else{
      numExps[nscStr] <- length(iNsc)
    }
  }
  
  return(numExps)
}