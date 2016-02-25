#' Returns a vector indicating the number of NCI-60 cell lines with missing activity data
#' for set of compounds.
#' 
#' @param nscSet a character vector specifying NSC identifier(s) for compound(s) of interest.
#' @return a numeric vector indicating the number of NCI-60 cell lines with missing 
#' activity data, indexed by the identifiers in nscSet.
#' 
#' @examples 
#' nscSet <- c("1", "17", "89", "609699")
#' getNumMissingLines(nscSet)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom Biobase exprs
getNumMissingLines <- function(nscSet){
	drugAct <- exprs(getAct(rcellminerData::drugData))
  
  nscSet <- as.character(nscSet)
  numMissingVals <- integer(length(nscSet))
  names(numMissingVals) <- nscSet
  
  for (nscStr in nscSet){
    numMissingVals[nscStr] <- sum(is.na(drugAct[nscStr, ]))
  }
  
  return(numMissingVals)
}